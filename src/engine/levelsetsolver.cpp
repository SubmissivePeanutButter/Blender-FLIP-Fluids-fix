/*
MIT License

Copyright (C) 2025 Ryan L. Guy & Dennis Fassbaender

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "levelsetsolver.h"

#include <limits>
#include <cmath>

#include "threadutils.h"
#include "grid3d.h"


LevelSetSolver::LevelSetSolver() {
}
struct Temp_Struct_4 {
    Array3d<float>* tempPtr;
    Array3d<float>* outputPtr;
    float dx;
    float dtau;
    std::vector<GridIndex>* solverCells;
    LevelSetSolver* Pointer;

};
void LevelSetSolver::reinitializeEno(Array3d<float> &inputSDF, 
                                     float dx, 
                                     float maxDistance, 
                                     std::vector<GridIndex> &solverCells, 
                                     Array3d<float> &outputSDF) {

    _maxCFL = 0.25;

    float dtau = _getPseudoTimeStep(inputSDF, dx);
    int numIterations = _getNumberOfIterations(maxDistance, dtau);

    int isize = inputSDF.width;
    int jsize = inputSDF.height;
    int ksize = inputSDF.depth;
    outputSDF = inputSDF;

    Array3d<float> tempSDF(isize, jsize, ksize);
    Array3d<float> *tempPtr = &tempSDF;
    Array3d<float> *outputPtr = &outputSDF;

    for (int n = 0; n < numIterations; n++) {

        int numCPU = ThreadUtils::getMaxThreadCount();
        int numthreads = (int)fmin(numCPU, solverCells.size());
        //std::vector<std::thread> threads(numthreads);
        //std::vector<int> intervals = ThreadUtils::splitRangeIntoIntervals(0, solverCells.size(), numthreads);
        //for (int i = 0; i < numthreads; i++) {
        //    threads[i] = std::thread(&LevelSetSolver::_stepSolverThreadEno, this,
        //                             intervals[i], intervals[i + 1], tempPtr, outputPtr, dx, dtau, &solverCells);
        //}

        //for (int i = 0; i < numthreads; i++) {
        //    threads[i].join();
        //}
        Temp_Struct_4 Temp{
            tempPtr,
            outputPtr,
            dx,
            dtau,
            &solverCells,
            this
        };
        ThreadUtils::Thread_Pool.Run_Function(_stepSolverThreadedEno, 0, solverCells.size(), &Temp, numthreads);
        ThreadUtils::Thread_Pool.Sync();
        std::swap(tempPtr, outputPtr);
    }

    float *rawarraysrc = outputPtr->getRawArray();
    float *rawarraydst = outputSDF.getRawArray();
    int n = outputSDF.getNumElements();
    for (int i = 0; i < n; i++) {
        rawarraydst[i] = rawarraysrc[i];
    }
}

void LevelSetSolver::reinitializeUpwind(Array3d<float> &inputSDF, 
                                        float dx, 
                                        float maxDistance, 
                                        std::vector<GridIndex> &solverCells, 
                                        Array3d<float> &outputSDF) {

    _maxCFL = 0.5;

    float dtau = _getPseudoTimeStep(inputSDF, dx);
    int numIterations = _getNumberOfIterations(maxDistance, dtau);

    int isize = inputSDF.width;
    int jsize = inputSDF.height;
    int ksize = inputSDF.depth;
    outputSDF = inputSDF;

    Array3d<float> tempSDF(isize, jsize, ksize);
    Array3d<float> *tempPtr = &tempSDF;
    Array3d<float> *outputPtr = &outputSDF;

    float lastMaxDiff = -1.0f;
    for (int n = 0; n < numIterations; n++) {
        int numCPU = ThreadUtils::getMaxThreadCount();
        int numthreads = (int)fmin(numCPU, solverCells.size());
        //std::vector<std::thread> threads(numthreads);
        //std::vector<int> intervals = ThreadUtils::splitRangeIntoIntervals(0, solverCells.size(), numthreads);
        //for (int i = 0; i < numthreads; i++) {
        //    threads[i] = std::thread(&LevelSetSolver::_stepSolverThreadUpwind, this,
        //                             intervals[i], intervals[i + 1], tempPtr, outputPtr, dx, dtau, &solverCells);
        //}

        //for (int i = 0; i < numthreads; i++) {
        //    threads[i].join();
        //}
        Temp_Struct_4 Temp{
            tempPtr,
            outputPtr,
            dx,
            dtau,
            &solverCells,
            this
        };
        ThreadUtils::Thread_Pool.Run_Function(LevelSetSolver::_stepSolverThreadedUpwind, 0, solverCells.size(), &Temp, numthreads);
        ThreadUtils::Thread_Pool.Sync();

        float maxDiff = 0;
        for (size_t cidx = 0; cidx < solverCells.size(); cidx++) {
            GridIndex g = solverCells[cidx];
            float diff = std::abs(tempPtr->get(g) - outputPtr->get(g));
            maxDiff = std::max(diff, maxDiff);
        }

        std::swap(tempPtr, outputPtr);

        if (std::abs(maxDiff - lastMaxDiff) < _upwindErrorThreshold * dx) {
            break;
        }
        lastMaxDiff = maxDiff;
    }

    float *rawarraysrc = outputPtr->getRawArray();
    float *rawarraydst = outputSDF.getRawArray();
    int n = outputSDF.getNumElements();
    for (int i = 0; i < n; i++) {
        rawarraydst[i] = rawarraysrc[i];
    }
}

float LevelSetSolver::_getPseudoTimeStep(Array3d<float> &sdf, float dx) {
    int isize = sdf.width;
    int jsize = sdf.height;
    int ksize = sdf.depth;

    float maxS = -std::numeric_limits<float>::max();
    float dtau = _maxCFL * dx;

    for (int k = 0; k < ksize; k++) {
        for (int j = 0; j < jsize; j++) {
            for (int i = 0; i < isize; i++) {
                float s = _sign(sdf, dx, i, j, k);
                maxS = std::max(s, maxS);
            }
        }
    }

    while (dtau * maxS / dx > _maxCFL) {
        dtau *= 0.5;
    }

    return dtau;
}

float LevelSetSolver::_sign(Array3d<float> &sdf, float dx, int i, int j, int k) {
    double d = sdf(i, j, k);
    return d / std::sqrt(d * d + dx * dx);
}

int LevelSetSolver::_getNumberOfIterations(float maxDistance, float dtau) {
    return static_cast<int>(std::ceil(maxDistance / dtau));
}

void LevelSetSolver::_stepSolverThreadEno(int startidx, int endidx, 
                                       Array3d<float> *tempPtr,
                                       Array3d<float> *outputPtr, 
                                       float dx,
                                       float dtau,
                                       std::vector<GridIndex> *solverCells) {

    std::array<float, 2> derx, dery, derz;
    for (int idx = startidx; idx < endidx; idx++) {
        GridIndex g = solverCells->at(idx);

        float s = _sign(*outputPtr, dx, g.i, g.j, g.k);
        _getDerivativesEno(outputPtr, g.i, g.j, g.k, dx, &derx, &dery, &derz);

        float val = outputPtr->get(g)
            - dtau * std::max(s, 0.0f)
                * (std::sqrt(_square(std::max(derx[0], 0.0f))
                           + _square(std::min(derx[1], 0.0f))
                           + _square(std::max(dery[0], 0.0f))
                           + _square(std::min(dery[1], 0.0f))
                           + _square(std::max(derz[0], 0.0f))
                           + _square(std::min(derz[1], 0.0f))) - 1.0f)
            - dtau * std::min(s, 0.0f)
                * (std::sqrt(_square(std::min(derx[0], 0.0f))
                           + _square(std::max(derx[1], 0.0f))
                           + _square(std::min(dery[0], 0.0f))
                           + _square(std::max(dery[1], 0.0f))
                           + _square(std::min(derz[0], 0.0f))
                           + _square(std::max(derz[1], 0.0f))) - 1.0f);

        tempPtr->set(g, val);
    }
}
void LevelSetSolver::_stepSolverThreadedEno(int startidx, int endidx, void* Data, int Thread_Number) {
    Temp_Struct_4* Temp = static_cast<Temp_Struct_4*>(Data);
    std::array<float, 2> derx, dery, derz;
    for (int idx = startidx; idx < endidx; idx++) {
        GridIndex g = Temp->solverCells->at(idx);

        float s = Temp->Pointer->_sign(*(Temp->outputPtr), Temp->dx, g.i, g.j, g.k);
        Temp->Pointer->_getDerivativesEno(Temp->outputPtr, g.i, g.j, g.k, Temp->dx, &derx, &dery, &derz);

        float val = Temp->outputPtr->get(g)
            - Temp->dtau * std::max(s, 0.0f)
            * (std::sqrt(Temp->Pointer->_square(std::max(derx[0], 0.0f))
                + Temp->Pointer->_square(std::min(derx[1], 0.0f))
                + Temp->Pointer->_square(std::max(dery[0], 0.0f))
                + Temp->Pointer->_square(std::min(dery[1], 0.0f))
                + Temp->Pointer->_square(std::max(derz[0], 0.0f))
                + Temp->Pointer->_square(std::min(derz[1], 0.0f))) - 1.0f)
            - Temp->dtau * std::min(s, 0.0f)
            * (std::sqrt(Temp->Pointer->_square(std::min(derx[0], 0.0f))
                + Temp->Pointer->_square(std::max(derx[1], 0.0f))
                + Temp->Pointer->_square(std::min(dery[0], 0.0f))
                + Temp->Pointer->_square(std::max(dery[1], 0.0f))
                + Temp->Pointer->_square(std::min(derz[0], 0.0f))
                + Temp->Pointer->_square(std::max(derz[1], 0.0f))) - 1.0f);

        Temp->tempPtr->set(g, val);
    }
}

void LevelSetSolver::_getDerivativesEno(Array3d<float> *grid,
                                        int i, int j, int k, float dx, 
                                        std::array<float, 2> *derx,
                                        std::array<float, 2> *dery,
                                        std::array<float, 2> *derz) {
    float D0[7];
    int isize = grid->width;
    int jsize = grid->height;
    int ksize = grid->depth;

    int im3 = (i < 3) ? 0 : i - 3;
    int im2 = (i < 2) ? 0 : i - 2;
    int im1 = (i < 1) ? 0 : i - 1;
    int ip1 = std::min(i + 1, isize - 1);
    int ip2 = std::min(i + 2, isize - 1);
    int ip3 = std::min(i + 3, isize - 1);
    int jm3 = (j < 3) ? 0 : j - 3;
    int jm2 = (j < 2) ? 0 : j - 2;
    int jm1 = (j < 1) ? 0 : j - 1;
    int jp1 = std::min(j + 1, jsize - 1);
    int jp2 = std::min(j + 2, jsize - 1);
    int jp3 = std::min(j + 3, jsize - 1);
    int km3 = (k < 3) ? 0 : k - 3;
    int km2 = (k < 2) ? 0 : k - 2;
    int km1 = (k < 1) ? 0 : k - 1;
    int kp1 = std::min(k + 1, ksize - 1);
    int kp2 = std::min(k + 2, ksize - 1);
    int kp3 = std::min(k + 3, ksize - 1);

    // 3rd-order ENO differencing
    D0[0] = grid->get(im3, j, k);
    D0[1] = grid->get(im2, j, k);
    D0[2] = grid->get(im1, j, k);
    D0[3] = grid->get(i, j, k);
    D0[4] = grid->get(ip1, j, k);
    D0[5] = grid->get(ip2, j, k);
    D0[6] = grid->get(ip3, j, k);
    *derx = _eno3(D0, dx);

    D0[0] = grid->get(i, jm3, k);
    D0[1] = grid->get(i, jm2, k);
    D0[2] = grid->get(i, jm1, k);
    D0[3] = grid->get(i, j, k);
    D0[4] = grid->get(i, jp1, k);
    D0[5] = grid->get(i, jp2, k);
    D0[6] = grid->get(i, jp3, k);
    *dery = _eno3(D0, dx);

    D0[0] = grid->get(i, j, km3);
    D0[1] = grid->get(i, j, km2);
    D0[2] = grid->get(i, j, km1);
    D0[3] = grid->get(i, j, k);
    D0[4] = grid->get(i, j, kp1);
    D0[5] = grid->get(i, j, kp2);
    D0[6] = grid->get(i, j, kp3);
    *derz = _eno3(D0, dx);
}

std::array<float, 2> LevelSetSolver::_eno3(float *D0, float dx) {
    float invdx = 1.0f / dx;
    float hinvdx = invdx / 2.0f;
    float tinvdx = invdx / 3.0f;
    float D1[6], D2[5], D3[2];
    float dQ1, dQ2, dQ3;
    float c, cstar;
    int Kstar;
    std::array<float, 2> dfx;

    D1[0] = invdx * (D0[1] - D0[0]);
    D1[1] = invdx * (D0[2] - D0[1]);
    D1[2] = invdx * (D0[3] - D0[2]);
    D1[3] = invdx * (D0[4] - D0[3]);
    D1[4] = invdx * (D0[5] - D0[4]);
    D1[5] = invdx * (D0[6] - D0[5]);

    D2[0] = hinvdx * (D1[1] - D1[0]);
    D2[1] = hinvdx * (D1[2] - D1[1]);
    D2[2] = hinvdx * (D1[3] - D1[2]);
    D2[3] = hinvdx * (D1[4] - D1[3]);
    D2[4] = hinvdx * (D1[5] - D1[4]);

    for (int K = 0; K < 2; K++) {
        if (std::fabs(D2[K + 1]) < std::fabs(D2[K + 2])) {
            c = D2[K + 1];
            Kstar = K - 1;
            D3[0] = tinvdx * (D2[K + 1] - D2[K]);
            D3[1] = tinvdx * (D2[K + 2] - D2[K + 1]);
        } else {
            c = D2[K + 2];
            Kstar = K;
            D3[0] = tinvdx * (D2[K + 2] - D2[K + 1]);
            D3[1] = tinvdx * (D2[K + 3] - D2[K + 2]);
        }

        if (std::fabs(D3[0]) < std::fabs(D3[1])) {
            cstar = D3[0];
        } else {
            cstar = D3[1];
        }

        dQ1 = D1[K + 2];
        dQ2 = c * (2 * (1 - K) - 1) * dx;
        dQ3 = cstar * (3 * ((1 - Kstar) * (1 - Kstar)) - 6 * (1 - Kstar) + 2) * dx * dx;

        dfx[K] = dQ1 + dQ2 + dQ3;
    }

    return dfx;
}

void LevelSetSolver::_stepSolverThreadUpwind(int startidx, int endidx, 
                                       Array3d<float> *tempPtr,
                                       Array3d<float> *outputPtr, 
                                       float dx,
                                       float dtau,
                                       std::vector<GridIndex> *solverCells) {

    std::array<float, 2> derx, dery, derz;
    for (int idx = startidx; idx < endidx; idx++) {
        GridIndex g = solverCells->at(idx);

        float s = _sign(*outputPtr, dx, g.i, g.j, g.k);
        _getDerivativesUpwind(outputPtr, g.i, g.j, g.k, dx, &derx, &dery, &derz);

        float val = outputPtr->get(g)
            - dtau * std::max(s, 0.0f)
                * (std::sqrt(_square(std::max(derx[0], 0.0f))
                           + _square(std::min(derx[1], 0.0f))
                           + _square(std::max(dery[0], 0.0f))
                           + _square(std::min(dery[1], 0.0f))
                           + _square(std::max(derz[0], 0.0f))
                           + _square(std::min(derz[1], 0.0f))) - 1.0f)
            - dtau * std::min(s, 0.0f)
                * (std::sqrt(_square(std::min(derx[0], 0.0f))
                           + _square(std::max(derx[1], 0.0f))
                           + _square(std::min(dery[0], 0.0f))
                           + _square(std::max(dery[1], 0.0f))
                           + _square(std::min(derz[0], 0.0f))
                           + _square(std::max(derz[1], 0.0f))) - 1.0f);

        tempPtr->set(g, val);
    }
}
void LevelSetSolver::_stepSolverThreadedUpwind(int startidx, int endidx, void* Data, int Thread_Number) {
    Temp_Struct_4* Temp = static_cast<Temp_Struct_4*>(Data);
    Array3d<float>* tempPtr = Temp->tempPtr;
    Array3d<float>* outputPtr = Temp->outputPtr;
    //temp
    int x1, y1, z1, bx1, by1, bz1, x2, y2, z2, bx2, by2, bz2;
    tempPtr->getGridDimensions(&x1, &y1, &z1);
    outputPtr->getGridDimensions(&x2, &y2, &z2);

    //
    
    float dx = Temp->dx;
    float dtau = Temp->dtau;
    std::vector<GridIndex>* solverCells = Temp->solverCells;
    LevelSetSolver* Pointer = Temp->Pointer;
    std::array<float, 2> derx, dery, derz;
    for (int idx = startidx; idx < endidx; idx++) {
        GridIndex g = solverCells->at(idx);
        //temp
        tempPtr->getGridDimensions(&bx1, &by1, &bz1);
        outputPtr->getGridDimensions(&bx2, &by2, &bz2);
        if ((bx1 != x1) || (by1 != y1) || (bz1 != z1)) {
            std::cout << "Hell";
        }
        if ((bx2 != x2) || (by2 != y2) || (bz2 != z2)) {
            std::cout << "Hell2";
        }
        //
        float s = Pointer->_sign(*outputPtr, dx, g.i, g.j, g.k);
        //temp
        tempPtr->getGridDimensions(&bx1, &by1, &bz1);
        outputPtr->getGridDimensions(&bx2, &by2, &bz2);
        if ((bx1 != x1) || (by1 != y1) || (bz1 != z1)) {
            std::cout << "Hell";
        }
        if ((bx2 != x2) || (by2 != y2) || (bz2 != z2)) {
            std::cout << "Hell2";
        }
        //
        Pointer->_getDerivativesUpwind(outputPtr, g.i, g.j, g.k, dx, &derx, &dery, &derz);
        //temp
        tempPtr->getGridDimensions(&bx1, &by1, &bz1);
        outputPtr->getGridDimensions(&bx2, &by2, &bz2);
        if ((bx1 != x1) || (by1 != y1) || (bz1 != z1)) {
            std::cout << "Hell";
        }
        if ((bx2 != x2) || (by2 != y2) || (bz2 != z2)) {
            std::cout << "Hell2";
        }
        //
        float val = outputPtr->get(g)
            - dtau * std::max(s, 0.0f)
            * (std::sqrt(Pointer->_square(std::max(derx[0], 0.0f))
                + Pointer->_square(std::max(dery[0], 0.0f))
                + Pointer->_square(std::min(derx[1], 0.0f))
                + Pointer->_square(std::min(dery[1], 0.0f))
                + Pointer->_square(std::max(derz[0], 0.0f))
                + Pointer->_square(std::min(derz[1], 0.0f))) - 1.0f)
            - dtau * std::min(s, 0.0f)
            * (std::sqrt(Pointer->_square(std::min(derx[0], 0.0f))
                + Pointer->_square(std::max(derx[1], 0.0f))
                + Pointer->_square(std::min(dery[0], 0.0f))
                + Pointer->_square(std::max(dery[1], 0.0f))
                + Pointer->_square(std::min(derz[0], 0.0f))
                + Pointer->_square(std::max(derz[1], 0.0f))) - 1.0f);
        //temp
        tempPtr->getGridDimensions(&bx1, &by1, &bz1);
        outputPtr->getGridDimensions(&bx2, &by2, &bz2);
        if ((bx1 != x1)||(by1 != y1)||(bz1 != z1)) {
            std::cout << "Hell";
        }
        if ((bx2 != x2) || (by2 != y2) || (bz2 != z2)) {
            std::cout << "Hell2";
        }
        //
        tempPtr->set(g, val);
        
    }

}

void LevelSetSolver::_getDerivativesUpwind(Array3d<float> *grid,
                                        int i, int j, int k, float dx, 
                                        std::array<float, 2> *derx,
                                        std::array<float, 2> *dery,
                                        std::array<float, 2> *derz) {
    
    float D0[3];
    int isize = grid->width;
    int jsize = grid->height;
    int ksize = grid->depth;

    int im1 = (i < 1) ? 0 : i - 1;
    int ip1 = std::min(i + 1, isize - 1);
    int jm1 = (j < 1) ? 0 : j - 1;
    int jp1 = std::min(j + 1, jsize - 1);
    int km1 = (k < 1) ? 0 : k - 1;
    int kp1 = std::min(k + 1, ksize - 1);

    D0[0] = grid->get(im1, j, k);
    D0[1] = grid->get(i, j, k);
    D0[2] = grid->get(ip1, j, k);
    *derx = _upwind1(D0, dx);

    D0[0] = grid->get(i, jm1, k);
    D0[1] = grid->get(i, j, k);
    D0[2] = grid->get(i, jp1, k);
    *dery = _upwind1(D0, dx);

    D0[0] = grid->get(i, j, km1);
    D0[1] = grid->get(i, j, k);
    D0[2] = grid->get(i, j, kp1);
    *derz = _upwind1(D0, dx);
}

std::array<float, 2> LevelSetSolver::_upwind1(float *D0, float dx) {
    float invdx = 1.0f / dx;
    std::array<float, 2> dfx;
    dfx[0] = invdx * (D0[1] - D0[0]);
    dfx[1] = invdx * (D0[2] - D0[1]);
    return dfx;
}