#include <boost/regex.hpp>

int main()
{
  boost::regex rx("a*b");

  if (boost::regex_match("aab", rx))
    ;
  return 0;
}
