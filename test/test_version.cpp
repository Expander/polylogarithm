#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "version.hpp"
#include <iostream>

TEST_CASE("test_version_string")
{
   using polylogarithm::version_major;
   using polylogarithm::version_minor;
   using polylogarithm::version_patch;

   std::cout << "Polylogarithm version: " << version_major << "."
             << version_minor << "." << version_patch << std::endl;
   CHECK(true);
}
