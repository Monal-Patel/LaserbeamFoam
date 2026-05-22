#pragma once
#include <chrono>
#include <polyMesh.H>
#ifdef LASERFOAM_PROFILING
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#endif

namespace laserfoam {

class AbstractTimer {
protected:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_;
  static double elapsed(const TimePoint &start, const TimePoint &end) {
    auto r = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    return r.count() * 1e-9;
  }
  TimePoint now() { return std::chrono::high_resolution_clock::now(); }
  AbstractTimer() { start_ = std::chrono::high_resolution_clock::now(); }
};

/** A timer for C++ blocks */
class BasicTimer : public AbstractTimer {
  std::string label_;

public:
  BasicTimer(std::string label) : label_(label) {}
  void printElapsed() {
    Foam::Info << label_.c_str() << ": " << elapsed(start_, now()) << "s\n";
  }
};

/** A timer for OpenFOAM time iteration */
class Chrono : public AbstractTimer {
  TimePoint prev_;
  double cumulatedCells_; // num cell * seconds
  const Foam::polyMesh &mesh_;
  int numIter_ = 0;

public:
  Chrono(const Foam::polyMesh &mesh) : mesh_(mesh) { prev_ = start_; }
  void nextIter(int numPimpleIter) {
    numIter_++;
    auto numCells = mesh_.globalData().nTotalCells() * numPimpleIter;
    auto current = now();
    auto duration = elapsed(prev_, current);
    auto totalDuration = elapsed(start_, current);
    cumulatedCells_ += numCells * duration;
    auto averageNumCells = cumulatedCells_ / totalDuration;
    double timeByCells = (Foam::UPstream::nProcs() * duration * 1e6) / numCells;
    double timeByCellsTotal = (Foam::UPstream::nProcs() * totalDuration * 1e6) /
                              averageNumCells / numIter_;
    Foam::Info
        << "Clock time by cell, iteration and processor: previous iteration="
        << timeByCells << " µs, average=" << timeByCellsTotal << " µs\n";
    prev_ = current;
  }
};
#ifdef LASERFOAM_PROFILING

class CumulativeTimer : public AbstractTimer {
  struct Counter {
    double seconds = 0.0;
    int calls = 0;
  };

  std::vector<std::string> order_;
  std::map<std::string, Counter> counters_;

public:
  CumulativeTimer(std::initializer_list<std::string> names) : order_(names) {
    for (const auto &name : order_)
      counters_[name] = {};
  }

  TimePoint tick() { return now(); }

  void accumulate(const std::string &name, TimePoint start) {
    auto &c = counters_[name];
    c.seconds += elapsed(start, now());
    c.calls++;
  }

  void report() {
    double total = 0.0;
    for (const auto &name : order_)
      total += counters_[name].seconds;
    Foam::Info << "--- Equation timing [s] (cumulative over PIMPLE) ---\n";
    for (const auto &name : order_) {
      const auto &c = counters_[name];
      char buf[128];
      std::snprintf(buf, sizeof(buf), "  %-12s %8.4f s  %3d calls  %3d%%",
                    name.c_str(), c.seconds, c.calls,
                    total > 0 ? int(100.0 * c.seconds / total + 0.5) : 0);
      Foam::Info << buf << "\n";
    }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "  %-12s %8.4f s", "TOTAL", total);
    Foam::Info << buf << "\n";
  }
};

#endif // LASERFOAM_PROFILING

} // namespace laserfoam
