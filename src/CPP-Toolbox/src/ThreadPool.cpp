/*
 * ThreadPool.cpp
 *
 */

#include "ThreadPool.h"

#if __cplusplus > 201103L

ThreadPool::ThreadPool(unsigned int noThreads, std::vector<SimInterface*> &vec) : refVecSim(vec), busyThreads(0), isRunning(true) {
    for(unsigned int i=0; i<noThreads; i++) {
        vecThreads.emplace_back(&ThreadPool::threadProcess, this);
    }
}

ThreadPool::~ThreadPool() {
    isRunning = false;
    cvWaitForTask.notify_all();
    for(auto iter = vecThreads.begin(); iter!=vecThreads.end(); ++iter)
        iter->join();
}

void ThreadPool::threadProcess() {
    while(isRunning) {
        std::unique_lock<std::mutex> lkTasks(mTasks);

        while(isRunning && tasks.empty()) {
            // wait for new task and prevent spurious wake-up
            cvWaitForTask.wait(lkTasks);//, [this]{ return !(isRunning && tasks.empty()); }
        }
        // check if still running
        if(!isRunning)
            return;

        ++busyThreads;

        // get and simulate new task
        Task t = tasks.front();
        tasks.pop();
        lkTasks.unlock();

//        std::cout << "Simulating " << t.simIndex << std::endl;
//        int tmp = std::floor( (double)std::rand()/RAND_MAX*300.0 );
//        std::this_thread::sleep_for(std::chrono::milliseconds  ( tmp ));
//        std::cout << "Ready simulating " << t.simIndex << std::endl;
        if(t.simIndex<refVecSim.size())
            *(t.success) = refVecSim[t.simIndex]->Simulate(t.obs, t.traj, nullptr);
        else
            *(t.success) = false;

        // push results
        {
            std::lock_guard<std::mutex> lkResults(mResults);
            results.push(t);
        }
        --busyThreads;
        cvWaitForTasksFinished.notify_one();
    }
}

void ThreadPool::pushTask(unsigned int simIndex, unsigned int resIndex, bool *success, double *obs, double *traj) {
    // careful here: somehow scoping this lock and destroying it
    // before the task notification might lead to a deadlock
    // in stress tests under mingw32 7.2.0
    std::lock_guard<std::mutex> lkTasks(mTasks);
    tasks.emplace(simIndex, resIndex, success, obs, traj);
    cvWaitForTask.notify_one();
}

void ThreadPool::clearTasks() {
    std::lock_guard<std::mutex> lkTasks(mTasks);
    std::queue<Task> empty;
    tasks.swap(empty);
}

unsigned int ThreadPool::waitForOne() {
    std::unique_lock<std::mutex> lkTasks(mTasks);
    std::unique_lock<std::mutex> lkResults(mResults);
    if(results.empty()) {
        lkResults.unlock(); // unlock for now
        while(!tasks.empty() || busyThreads>0) {
            cvWaitForTasksFinished.wait(lkTasks);//, [this]{return !(!tasks.empty() || busyThreads>0);}
        }
        lkResults.lock(); // reacquire lock to pop result
    }
    if(results.empty()) { // still empty...no tasks and busy threads
        return UINT_MAX;
    }
    unsigned int index = results.front().resIndex;
    results.pop();
    return index;
}

void ThreadPool::waitForAll() {
    while(waitForOne()!=UINT_MAX) {}
//    // stress test code:
//    std::set<unsigned int> allRes;
//    while(true) {
//        unsigned int tmp = waitForOne();
//        if(tmp==UINT_MAX)
//            break;
//        else
//            allRes.insert(tmp);
//    }
//    std::cout << "Total: " << allRes.size() << std::endl;
//    for(unsigned int i=0; i<1000; i++) {
//        size_t no = allRes.erase(i);
//        if(no==0)
//            std::cout << "Missing: " << i << std::endl;
//    }
}

#endif
