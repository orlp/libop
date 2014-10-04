#define CATCH_CONFIG_MAIN
#include "catch.h"

#include "../op.h"


TEST_CASE("op::pow") {
    REQUIRE(op::pow<0>(2) == 1);
    REQUIRE(op::pow<1>(2) == 2);
    REQUIRE(op::pow<2>(2) == 4);
    REQUIRE(op::pow<3>(2) == 8);
    REQUIRE(op::pow<0>(32) == 1);
    REQUIRE(op::pow<32>(0) == 0);
    REQUIRE(op::pow<-1>(2.) == Approx(0.5));
    REQUIRE(op::pow<2>(-1) == 1);
}


TEST_CASE("op::ilog") {
    REQUIRE(op::ilog<2>(-1) == -1);
    REQUIRE(op::ilog<2>(0) == -1);
    REQUIRE(op::ilog<2>(1) == 0);
    REQUIRE(op::ilog<2>(2) == 1);
    REQUIRE(op::ilog<2>(3) == 1);
    REQUIRE(op::ilog<2>(4) == 2);
    REQUIRE(op::ilog<2>(32) == 5);
    REQUIRE(op::ilog<3>(2) == 0);
    REQUIRE(op::ilog<3>(3) == 1);
    REQUIRE(op::ilog<3>(4) == 1);
    REQUIRE(op::ilog<3>(8) == 1);
    REQUIRE(op::ilog<3>(9) == 2);
}

