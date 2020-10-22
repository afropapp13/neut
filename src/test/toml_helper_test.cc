#include "toml/toml_helper.h"

#include <array>
#include <iostream>
#include <string>

using namespace toml::literals::toml_literals;

int main() {
  toml::value const &tab = R"(
		string = "string"
    bool_str = "on"
    float = 1
    float_vector = [1,1.5,6.0]
    [table]
    my_str = "string"
    [table.nested]
    my_str = "string"
	)"_toml;

  std::cout << toml_h::find<std::string>(tab,"string") << std::endl;
  std::cout << toml_h::find<bool>(tab, "bool_str") << std::endl;
  std::cout << toml_h::find<float>(tab, "float") << std::endl;

  auto float_vector = toml_h::find<std::vector<float> >(tab, "float_vector");
  auto float_array = toml_h::find<std::array<float, 3> >(tab, "float_vector");

  auto table = toml_h::find(tab, "table");
  auto nested = toml_h::find(table, "nested");
  auto nested_direct = toml_h::find_rec(tab, {"table", "nested"});
  auto nested_direct_str =
      toml_h::find_rec<std::string>(tab, {"table", "nested", "my_str"});
}