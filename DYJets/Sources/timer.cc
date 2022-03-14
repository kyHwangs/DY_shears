#include "timer.h"

#include <cmath>
#include <iomanip>
#include <iostream>

namespace util
{

timer::timer(timer::counter last, timer::counter begin)
    : _steps(last), _from(begin), _done(begin), _stop(false)
{
}

timer::~timer()
{
    if (_running) {
        stop();
    }
}

namespace /* anonymous */
{

static void print_time(timer::counter us)
{
    using namespace std;
    using chrono::duration_cast;

    chrono::microseconds time(us);

    int hours = duration_cast<chrono::hours>(time).count();
    if (hours > 0) {
        std::cerr << hours << "h ";
        time -= chrono::hours(hours);
    }
    int minutes = duration_cast<chrono::minutes>(time).count();
    if (minutes > 0) {
        std::cerr << std::setw(2) << minutes << "min ";
        time -= chrono::minutes(minutes);
    }
    if (hours == 0 && minutes < 5) {
        // Use milliseconds to round correctly
        int milliseconds = duration_cast<chrono::milliseconds>(time).count();
        std::cerr << std::setw(2) << std::ceil(milliseconds / 1000.) << "s";
    }
}
} // namespace anonymous

void timer::displayprogress(const timer::time_point &now, bool at_end) const
{
    using us = std::chrono::microseconds;
    using std::chrono::duration_cast;

    // Clear line
    std::cerr << "\033[2K";

    // Percentage
    std::cerr << std::fixed << std::setprecision(0) << std::setw(3) << std::right
              << 100. * (_done - _from) / (_steps - _from) << "%\t";

    // done/total
    std::cerr << " " << std::setw(std::ceil(std::log10(steps()))) << _done << "/" << _steps << "\t";

    if (_done - _from > (steps() - _from) / 10 || now - _start > std::chrono::seconds(2)) {
        // Time computations (in microseconds)
        us::rep totalduration = duration_cast<us>(now - _start).count();
        us::rep stepduration = duration_cast<us>(now - _last).count();

        // Speeds (s^-1)
        double instspeed = 1e6 * (_done - _lastcounter) / stepduration;
        double avgspeed = 1e6 * (_done - _from) / totalduration;

        if (!at_end) {
            std::cerr << " Speed: " << instspeed << "/s"
                      << " (average " << avgspeed << "/s)\t";
        } else {
            std::cerr << " Average speed: " << avgspeed << "/s"
                      << " (" << 1e6 / avgspeed << "us)\t";
        }

        // ETA
        if (!at_end) {
            std::cerr << " Time left: ";
            print_time(1e6 * (steps() - _from) / avgspeed - totalduration);
        } else {
            std::cerr << " Duration: ";
            print_time(totalduration);
        }
    }

    std::cerr << std::setprecision(6);
    std::cerr << "\r" << std::flush; // Back to beginning of line
}

void timer::start()
{
    if (_running) {
        return;
    }
    _running = true;
    _thread = std::thread([&]() {
        // Main work function, will be invoked in a separate thread
        _start = clock::now();
        _last = _start;

        displayprogress(_start);

        do { // Thread main loop
            std::this_thread::sleep_for(std::chrono::seconds(1));

            time_point now = clock::now();
            displayprogress(now);

            _last = now;
            _lastcounter = _done;
        } while (!_stop);

        time_point now = clock::now();
        displayprogress(now, true);

        std::cerr << std::endl;
    });
}

void timer::stop()
{
    _stop = true;
    _thread.join();
    _running = false;
}
} // namespace util
