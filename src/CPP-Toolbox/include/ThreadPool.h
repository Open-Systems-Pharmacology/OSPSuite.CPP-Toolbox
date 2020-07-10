/*
 * SimThreads.h
 *
 *  Created on: 15.09.2016
 *      Author: GGKMT
 */

#ifndef SRC_THREADPOOL_H_
#define SRC_THREADPOOL_H_

#if __cplusplus < 201103L
#warning C++11 support not available, disabling threads
#else

#include <vector>
#include <queue>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>

#include "Integrator.h"

struct Task {
    unsigned int simIndex = 0;
    unsigned int resIndex = 0;
    bool *success = nullptr;
    double *obs   = nullptr;
    double *traj  = nullptr;

    Task(unsigned int _simIndex, unsigned int _resIndex, bool *_success, double *_obs, double *_traj) :
        simIndex(_simIndex), resIndex(_resIndex), success(_success), obs(_obs), traj(_traj) {}
};

class ThreadPool {
private:
    std::vector<std::thread> vecThreads;

    std::queue<Task> tasks;
    std::queue<Task> results;

    std::mutex mTasks;
    std::mutex mResults;

    std::condition_variable cvWaitForTask;
    std::condition_variable cvWaitForTasksFinished;

    const std::vector<SimInterface*> &refVecSim;

    // determine if any threads are still busy (corresponding task
    // would be not in tasks anymore, but not in results yet)
    std::atomic<unsigned int> busyThreads;
    std::atomic<bool> isRunning;

    //	void init(unsigned int noThreads);
    void threadProcess();

public:
    //	SimThreadPool(unsigned int noThreads, SimInterface* s);
    ThreadPool(unsigned int noThreads, std::vector<SimInterface*> &vec);
    ~ThreadPool();
    ThreadPool(const ThreadPool&)            = delete;
    ThreadPool(ThreadPool&&)                 = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&)      = delete;

    void pushTask(unsigned int simIndex, unsigned int resIndex, bool *success, double *obs, double *traj);
    void clearTasks();
    unsigned int waitForOne();
    void waitForAll();
};

#endif /* __cplusplus version */
#endif /* SRC_THREADPOOL_H_ */
