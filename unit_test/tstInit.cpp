
#include <Kokkos_Core.hpp>

#include "CAinitialize.hpp"

#include <gtest/gtest.h>

#include "mpi.h"

#include <fstream>
#include <string>
#include <vector>

namespace Test {
//---------------------------------------------------------------------------//
void testReadTemperatureData() {}

//---------------------------------------------------------------------------//
// RUN TESTS
//---------------------------------------------------------------------------//
TEST(TEST_CATEGORY, temperature_init_test) { testReadTemperatureData(); }

} // end namespace Test
