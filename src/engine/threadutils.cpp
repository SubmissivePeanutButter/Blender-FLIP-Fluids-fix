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

#include <iostream>

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
	Number_Of_Threads_For_Current_Task = Number_Of_Threads;
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
void ThreadUtils::Thread_Pool_Handeler::Run_Function(std::function<void(int Start_Index, int End_Index, void* Data, int Thread_Number)> Task, int Range_Start, int Range_End, void* Data, int Number_Of_Threads_For_Task) {
	//handle stupidity 
	if (Number_Of_Threads_For_Task > 0) {// }:(
		//wait for threads finish
		for (unsigned int i = 0; i < Number_Of_Threads; i++) {
			Thread_Finished_Flag_Array[i].wait(false);
			Thread_Finished_Flag_Array[i].clear();
		}
		//updata how mant threads to use for task
		Number_Of_Threads_For_Current_Task = Number_Of_Threads_For_Task;

		//update task
		Thread_Task = Task;

		//update task data
		Task_Data = Data;

		//update task range
		Task_Range_Start = Range_Start;
		Task_Range_End = Range_End;

		//flip the flip flag
		if (Flip_To_Notify_Threads_Flag.test()) {
			Flip_To_Notify_Threads_Flag.clear();
		}
		else {
			Flip_To_Notify_Threads_Flag.test_and_set();
		}
		Flip_To_Notify_Threads_Flag.notify_all();
	}
	else{
		Sync();
	}
	return;
}
void ThreadUtils::Thread_Pool_Handeler::Dummy_Function(int Start_Index, int End_Index, void* Data, int Thread_Number) {
	return;
}
void ThreadUtils::Thread_Pool_Handeler::Sync() {
	//wait for threads finish
	for (unsigned int i = 0; i < Number_Of_Threads; i++) {
		Thread_Finished_Flag_Array[i].wait(false);
		//Thread_Finished_Flag_Array[i].test_and_set();
	}
	return;
}
	
void ThreadUtils::Thread_Pool_Handeler::Thread_Manager_Function(ThreadUtils::Thread_Pool_Handeler* Thread_Pool_Handler, unsigned int Thread_Number) {
	bool Flag_Proxy = Thread_Pool_Handler->Flip_To_Notify_Threads_Flag.test();
	//Thread_Number has range from 0 to number of threads -1
	do {
		//calculate task range
		int Total_Tasks_To_Run = Thread_Pool_Handler->Task_Range_End - Thread_Pool_Handler->Task_Range_Start;

		//same as interval size
		//int Number_Of_Tasks_Per_Thread = std::floor(Total_Tasks_To_Run / Thread_Pool_Handler->Number_Of_Threads);
		int Number_Of_Tasks_Per_Thread = std::floor(Total_Tasks_To_Run / Thread_Pool_Handler->Number_Of_Threads_For_Current_Task);// if 0, :(

		//same as interval remainder
		//int Task_Remainder = Total_Tasks_To_Run - (Number_Of_Tasks_Per_Thread * Thread_Pool_Handler->Number_Of_Threads);
		int Task_Remainder = Total_Tasks_To_Run - (Number_Of_Tasks_Per_Thread * Thread_Pool_Handler->Number_Of_Threads_For_Current_Task);

		int Thread_Start_Index = Thread_Pool_Handler->Task_Range_Start + (Number_Of_Tasks_Per_Thread * Thread_Number);

		int Thread_End_Index = Thread_Start_Index + Number_Of_Tasks_Per_Thread - 1;

		//distribute remaining tasks across threads
		if (Task_Remainder > 0) {//if less than 0, hell froze over
			//shift tasks over
			//use std::min to only shift by added space
			Thread_Start_Index += std::min(Thread_Number, (unsigned int)Task_Remainder);
			Thread_End_Index += std::min(Thread_Number, (unsigned int)Task_Remainder);
			//shift in space for remaining tasks tasks
			if (Thread_Number < Task_Remainder) {
				//handle all tasks being remainders
				if (Number_Of_Tasks_Per_Thread > 0) {//if less than 0, hell froze over again
					Thread_End_Index++;
				}
			}
		}
		//debug
		////if (Thread_Number == (Thread_Pool_Handler->Number_Of_Threads - 1)) {
		//if (Thread_Number == (Thread_Pool_Handler->Number_Of_Threads_For_Current_Task - 1)) {
		//	std::cout << "Total_Tasks_To_Run = " << Total_Tasks_To_Run << " Number_Of_Tasks_Per_Thread = " << Number_Of_Tasks_Per_Thread << " Task_Difference = " << Task_Difference << " Thread_Start_Index = " << Thread_Start_Index << " Thread_End_Index = " << Thread_End_Index << "\n";
		//}
		
		//handle having no tasks
		if (Total_Tasks_To_Run > 0) {//if less than 0, hell should turn into a ski resort
			//handle having more threads than threads per task
			if (Thread_Number < Thread_Pool_Handler->Number_Of_Threads_For_Current_Task) {
				//execute task
				Thread_Pool_Handler->Thread_Task(Thread_Start_Index, Thread_End_Index, Thread_Pool_Handler->Task_Data, Thread_Number);
			}
		}


		//set flag when task is finished
		Thread_Pool_Handler->Thread_Finished_Flag_Array[Thread_Number].test_and_set();
		Thread_Pool_Handler->Thread_Finished_Flag_Array[Thread_Number].notify_all();

		//wait for the flip flag
		Thread_Pool_Handler->Flip_To_Notify_Threads_Flag.wait(Flag_Proxy);
		Flag_Proxy = !Flag_Proxy;

		//exit whed done
	} while (!Thread_Pool_Handler->Exit_Flag.test());
	return;
}