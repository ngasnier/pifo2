#define BOOST_TEST_DYN_LINK  
#define BOOST_TEST_MODULE PifoTestcases

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(PifoTest) {
  BOOST_CHECK_EQUAL(1., 1.);
}
