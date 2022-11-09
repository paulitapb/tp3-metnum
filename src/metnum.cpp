#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include "eigen.h"

using namespace std;

namespace py=pybind11;

PYBIND11_MODULE(metnum, m) {
    m.def(
        "power_iteration", &power_iteration,
        "doingthings"
    );
    m.def(
        "deflation", &deflation,
        "doingthings"
    );
}