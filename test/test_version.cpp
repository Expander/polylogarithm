#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "version.h"
#include <iostream>

TEST_CASE("test_version_string")
{
   std::cout << "Polylogarithm version: " << POLYLOGARITHM_MAJOR << "."
             << POLYLOGARITHM_MINOR << "." << POLYLOGARITHM_PATCH << std::endl;
   CHECK(true);
}
