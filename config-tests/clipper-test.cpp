// test to see if the clipper version is up-to-date enough

#include <clipper/core/atomsf.h>

int main()
{
  auto sf = clipper::ScatteringFactors::instance()["C"];
  return 0;
}
