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

#include "threadutils.h"

#include <cmath>

#include "fluidsimassert.h"

int ThreadUtils::_maxThreadCount = 0;
bool ThreadUtils::_isMaxThreadCountInitialized = false;
ThreadUtils::Thread_Pool_Handeler ThreadUtils::Thread_Pool;
void ThreadUtils::_initializeMaxThreadCount() {
	if (_isMaxThreadCountInitialized) {
		return;
	}

	_maxThreadCount = std::thread::hardware_concurrency();
	_isMaxThreadCountInitialized = true;
}

int ThreadUtils::getMaxThreadCount(){
	_initializeMaxThreadCount();
	return _maxThreadCount;
}

void ThreadUtils::setMaxThreadCount(int n) {
	FLUIDSIM_ASSERT(n > 0);
	_maxThreadCount = n;
	_isMaxThreadCountInitialized = true;
}

std::vector<int> ThreadUtils::splitRangeIntoIntervals(int rangeBegin, int rangeEnd, 
													  int numIntervals) {
	int intervalSize = floor((double)(rangeEnd - rangeBegin) / (double)numIntervals);
	int intervalRemainder = (rangeEnd - rangeBegin) - intervalSize * numIntervals;
	std::vector<int> intervals;
	intervals.reserve(numIntervals + 1);
	intervals.push_back(rangeBegin);

	int intervalBegin = rangeBegin;
	for (int i = 0; i < numIntervals; i++) {
		int intervalEnd = intervalBegin + intervalSize;
		if (i < intervalRemainder) {
			intervalEnd++;
		}
		intervals.push_back(intervalEnd);
		intervalBegin = intervalEnd;
	}

	return intervals;
}
ThreadUtils::Thread_Pool_Handeler::Thread_Pool_Handeler() {
	Number_Of_Threads = std::thread::hardware_concurrency();
	Flip_To_Notify_Threads_Flag.clear();
	Thread_Finished_Flag_Array = new std::atomic_flag[Number_Of_Threads];
	Thread_Task = Dummy_Function;
	Task_Data = nullptr;
	Task_Range_Start = 0;
	Task_Range_End = 0;
	Thread_Array = static_cast<std::thread*>( operator new(sizeof(std::thread) * Number_Of_Threads));
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Finished_Flag_Array[i].clear();
		new (&Thread_Array[i])std::thread(Thread_Manager_Function, this, i);
	}
}
ThreadUtils::Thread_Pool_Handeler::~Thread_Pool_Handeler() {
	//wait for threads finish
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Finished_Flag_Array[i].wait(false);
	}
	//notify threads to close
	Exit_Flag.test_and_set();
	Exit_Flag.notify_all();

	//flip the flip flag
	if (Flip_To_Notify_Threads_Flag.test()) {
		Flip_To_Notify_Threads_Flag.clear();
	}
	else {
		Flip_To_Notify_Threads_Flag.test_and_set();
	}
	Flip_To_Notify_Threads_Flag.notify_all();

	//join threads
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Array[i].join();
	}

	//delete Thread_Array
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Array[i].~thread();
	}
	delete Thread_Array;

	//delete Thread_Finished_Flag_Array
	delete[] Thread_Finished_Flag_Array;
}
void ThreadUtils::Thread_Pool_Handeler::Run_Function(std::function<void(int Start_Index, int End_Index, void* Data)> Task, int Range_Start, int Range_End, void* Data) {
	//wait for threads finish
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Finished_Flag_Array[i].wait(false);
		Thread_Finished_Flag_Array[i].clear();
	}
	//update task
	Thread_Task = Task;

	//update task data
	Task_Data = Data;

	//flip the flip flag
	if (Flip_To_Notify_Threads_Flag.test()) {
		Flip_To_Notify_Threads_Flag.clear();
	}
	else {
		Flip_To_Notify_Threads_Flag.test_and_set();
	}
	Flip_To_Notify_Threads_Flag.notify_all();
}
void ThreadUtils::Thread_Pool_Handeler::Dummy_Function(int Start_Index, int End_Index, void* Data) {
	return;
}
void ThreadUtils::Thread_Pool_Handeler::Thread_Manager_Function(ThreadUtils::Thread_Pool_Handeler* Thread_Pool_Handler, unsigned int i) {
	bool Flag_Proxy = Thread_Pool_Handler->Flip_To_Notify_Threads_Flag.test();
	int Total_Tasks_To_Run = Thread_Pool_Handler->Task_Range_End - Thread_Pool_Handler->Task_Range_Start;
	int Number_Of_Tasks_Per_Thread = std::floor(Total_Tasks_To_Run / Thread_Pool_Handler->Number_Of_Threads);
	int Task_Difference = Total_Tasks_To_Run - (Number_Of_Tasks_Per_Thread * Thread_Pool_Handler->Number_Of_Threads);
	int Thread_Start_Index = Number_Of_Tasks_Per_Thread * (Thread_Pool_Handler->Number_Of_Threads - 1);
	int Thread_End_Index = (Thread_Start_Index + Number_Of_Tasks_Per_Thread) - 1;

	if ( i == Thread_Pool_Handler->Number_Of_Threads) {
		Thread_End_Index++;
		Thread_End_Index += Task_Difference;
	}
	do {
		//execute task
		Thread_Pool_Handler->Thread_Task(Thread_Start_Index, Thread_End_Index, Thread_Pool_Handler->Task_Data);

		//set flag when task is finished
		Thread_Pool_Handler->Thread_Finished_Flag_Array[i].test_and_set();
		Thread_Pool_Handler->Thread_Finished_Flag_Array[i].notify_all();

		//wait for the flip flag
		Thread_Pool_Handler->Flip_To_Notify_Threads_Flag.wait(Flag_Proxy);
		Flag_Proxy = !Flag_Proxy;

		//exit whed done
	} while (!Thread_Pool_Handler->Exit_Flag.test());
	return;
}