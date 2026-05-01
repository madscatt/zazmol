#pragma once

#include "sasmol/types.hpp"

#include <map>
#include <string>

namespace sasmol {

[[nodiscard]] const std::map<std::string, calc_type>& amu();

}  // namespace sasmol
