#include "sasmol/selection.hpp"

#include <cctype>
#include <charconv>
#include <functional>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace sasmol {
namespace {

enum class TokenKind {
  identifier,
  integer,
  string_literal,
  left_paren,
  right_paren,
  left_bracket,
  right_bracket,
  comparison,
  logical_and,
  logical_or,
  logical_not,
  null_literal,
  bool_literal,
  end,
};

struct Token {
  TokenKind kind{TokenKind::end};
  std::string text;
};

struct Predicate {
  std::function<bool(std::size_t)> eval;
};

std::string selection_error(const std::string& message) {
  return "selection parse failed: " + message;
}

std::string lower_ascii(std::string value) {
  for (auto& character : value) {
    character = static_cast<char>(
        std::tolower(static_cast<unsigned char>(character)));
  }
  return value;
}

std::string trim(std::string_view value) {
  auto begin = value.begin();
  auto end = value.end();
  while (begin != end && std::isspace(static_cast<unsigned char>(*begin))) {
    ++begin;
  }
  while (begin != end && std::isspace(static_cast<unsigned char>(*(end - 1)))) {
    --end;
  }
  return {begin, end};
}

bool is_identifier_start(char value) {
  return std::isalpha(static_cast<unsigned char>(value)) || value == '_';
}

bool is_identifier_char(char value) {
  return std::isalnum(static_cast<unsigned char>(value)) || value == '_';
}

std::optional<int> parse_int(std::string_view text) {
  int value{};
  const auto* begin = text.data();
  const auto* end = text.data() + text.size();
  const auto result = std::from_chars(begin, end, value);
  if (result.ec != std::errc{} || result.ptr != end) {
    return std::nullopt;
  }
  return value;
}

std::vector<Token> tokenize(const std::string& expression) {
  std::vector<Token> tokens;
  for (std::size_t offset = 0; offset < expression.size();) {
    const auto value = expression[offset];
    if (std::isspace(static_cast<unsigned char>(value))) {
      ++offset;
      continue;
    }
    if (is_identifier_start(value)) {
      const auto start = offset++;
      while (offset < expression.size() && is_identifier_char(expression[offset])) {
        ++offset;
      }
      auto text = expression.substr(start, offset - start);
      if (text == "and") {
        tokens.push_back({TokenKind::logical_and, text});
      } else if (text == "or") {
        tokens.push_back({TokenKind::logical_or, text});
      } else if (text == "not") {
        tokens.push_back({TokenKind::logical_not, text});
      } else if (text == "None") {
        tokens.push_back({TokenKind::null_literal, text});
      } else if (text == "True" || text == "False") {
        tokens.push_back({TokenKind::bool_literal, text});
      } else {
        tokens.push_back({TokenKind::identifier, text});
      }
      continue;
    }
    if (std::isdigit(static_cast<unsigned char>(value)) || value == '-') {
      const auto start = offset++;
      while (offset < expression.size() &&
             std::isdigit(static_cast<unsigned char>(expression[offset]))) {
        ++offset;
      }
      tokens.push_back({TokenKind::integer, expression.substr(start, offset - start)});
      continue;
    }
    if (value == '"' || value == '\'') {
      const auto quote = value;
      const auto start = ++offset;
      while (offset < expression.size() && expression[offset] != quote) {
        ++offset;
      }
      if (offset == expression.size()) {
        throw std::invalid_argument(selection_error("unterminated string literal"));
      }
      tokens.push_back(
          {TokenKind::string_literal, expression.substr(start, offset - start)});
      ++offset;
      continue;
    }
    if (value == '(') {
      tokens.push_back({TokenKind::left_paren, "("});
      ++offset;
      continue;
    }
    if (value == ')') {
      tokens.push_back({TokenKind::right_paren, ")"});
      ++offset;
      continue;
    }
    if (value == '[') {
      tokens.push_back({TokenKind::left_bracket, "["});
      ++offset;
      continue;
    }
    if (value == ']') {
      tokens.push_back({TokenKind::right_bracket, "]"});
      ++offset;
      continue;
    }
    if (value == '=' || value == '!' || value == '<' || value == '>') {
      const auto start = offset++;
      if (offset < expression.size() && expression[offset] == '=') {
        ++offset;
      }
      auto text = expression.substr(start, offset - start);
      if (text == "=") {
        text = "==";
      }
      if (text != "==" && text != "!=" && text != "<" && text != "<=" &&
          text != ">" && text != ">=") {
        throw std::invalid_argument(selection_error("unsupported comparison '" +
                                                    text + "'"));
      }
      tokens.push_back({TokenKind::comparison, text});
      continue;
    }
    throw std::invalid_argument(selection_error("unsupported token near '" +
                                                expression.substr(offset, 1) + "'"));
  }
  tokens.push_back({TokenKind::end, ""});
  return tokens;
}

class Parser {
 public:
  Parser(const Molecule& molecule, std::vector<Token> tokens)
      : molecule_(molecule), tokens_(std::move(tokens)) {}

  Predicate parse() {
    auto predicate = parse_or();
    expect(TokenKind::end, "end of expression");
    return predicate;
  }

 private:
  [[nodiscard]] const Token& current() const { return tokens_[position_]; }

  bool match(TokenKind kind) {
    if (current().kind != kind) {
      return false;
    }
    ++position_;
    return true;
  }

  Token expect(TokenKind kind, const char* what) {
    if (current().kind != kind) {
      throw std::invalid_argument(selection_error(std::string("expected ") + what));
    }
    return tokens_[position_++];
  }

  Predicate parse_or() {
    auto left = parse_and();
    while (match(TokenKind::logical_or)) {
      auto right = parse_and();
      auto left_eval = std::move(left.eval);
      auto right_eval = std::move(right.eval);
      left.eval = [left_eval, right_eval](std::size_t atom) {
        return left_eval(atom) || right_eval(atom);
      };
    }
    return left;
  }

  Predicate parse_and() {
    auto left = parse_unary();
    while (match(TokenKind::logical_and)) {
      auto right = parse_unary();
      auto left_eval = std::move(left.eval);
      auto right_eval = std::move(right.eval);
      left.eval = [left_eval, right_eval](std::size_t atom) {
        return left_eval(atom) && right_eval(atom);
      };
    }
    return left;
  }

  Predicate parse_unary() {
    if (match(TokenKind::logical_not)) {
      auto predicate = parse_unary();
      auto eval = std::move(predicate.eval);
      return {[eval](std::size_t atom) { return !eval(atom); }};
    }
    return parse_primary();
  }

  Predicate parse_primary() {
    if (match(TokenKind::left_paren)) {
      auto predicate = parse_or();
      expect(TokenKind::right_paren, "')'");
      return predicate;
    }
    return parse_comparison();
  }

  Predicate parse_comparison() {
    const auto field = expect(TokenKind::identifier, "descriptor name").text;
    expect(TokenKind::left_bracket, "'['");
    const auto index_name = expect(TokenKind::identifier, "atom index variable").text;
    if (index_name != "i") {
      throw std::invalid_argument(
          selection_error("only atom index variable 'i' is supported"));
    }
    expect(TokenKind::right_bracket, "']'");
    std::optional<int> string_offset;
    if (match(TokenKind::left_bracket)) {
      const auto offset = expect(TokenKind::integer, "string offset").text;
      const auto parsed_offset = parse_int(offset);
      if (!parsed_offset || *parsed_offset < 0) {
        throw std::invalid_argument(selection_error("invalid string offset"));
      }
      string_offset = *parsed_offset;
      expect(TokenKind::right_bracket, "']'");
    }
    const auto op = expect(TokenKind::comparison, "comparison operator").text;
    const auto value = current();
    if (value.kind != TokenKind::string_literal &&
        value.kind != TokenKind::integer && value.kind != TokenKind::null_literal &&
        value.kind != TokenKind::bool_literal) {
      throw std::invalid_argument(
          selection_error("comparison value must be a string, integer, None, or bool literal"));
    }
    ++position_;

    if (field == "record" || field == "name" || field == "loc" ||
        field == "resname" || field == "chain" || field == "rescode" ||
        field == "occupancy" || field == "beta" || field == "segname" ||
        field == "element" || field == "charge" || field == "moltype" ||
        field == "charmm_type") {
      if (value.kind != TokenKind::string_literal &&
          value.kind != TokenKind::null_literal) {
        throw std::invalid_argument(
            selection_error("descriptor '" + field + "' requires a string or None literal"));
      }
      return string_comparison(field, op, value.text, string_offset,
                               value.kind == TokenKind::null_literal);
    }
    if (string_offset) {
      throw std::invalid_argument(
          selection_error("string offsets are only supported for string descriptors"));
    }
    if (field == "resid" || field == "index" || field == "original_index" ||
        field == "original_resid" || field == "residue_flag") {
      if (value.kind != TokenKind::integer && value.kind != TokenKind::bool_literal) {
        throw std::invalid_argument(
            selection_error("descriptor '" + field + "' requires an integer or bool literal"));
      }
      int expected{};
      if (value.kind == TokenKind::bool_literal) {
        expected = value.text == "True" ? 1 : 0;
      } else {
        const auto parsed = parse_int(value.text);
        if (!parsed) {
          throw std::invalid_argument(selection_error("invalid integer literal"));
        }
        expected = *parsed;
      }
      return int_comparison(field, op, expected);
    }
    throw std::invalid_argument(selection_error("unsupported descriptor '" + field + "'"));
  }

  Predicate string_comparison(const std::string& field, const std::string& op,
                              const std::string& expected,
                              std::optional<int> string_offset,
                              bool expected_none) const {
    if (op != "==" && op != "!=") {
      if (expected_none) {
        throw std::invalid_argument(
            selection_error("None comparisons only support == and !="));
      }
    }
    const auto& values = string_descriptor(field);
    require_length(field, values.size());
    return {[values = &values, op, expected, string_offset,
             expected_none](std::size_t atom) {
      bool equal = false;
      if (!expected_none) {
        auto actual = (*values)[atom];
        if (string_offset) {
          const auto offset = static_cast<std::size_t>(*string_offset);
          actual = offset < actual.size() ? actual.substr(offset, 1) : std::string{};
        }
        equal = actual == expected;
        if (op == "<") return actual < expected;
        if (op == "<=") return actual <= expected;
        if (op == ">") return actual > expected;
        if (op == ">=") return actual >= expected;
      }
      return op == "==" ? equal : !equal;
    }};
  }

  Predicate int_comparison(const std::string& field, const std::string& op,
                           int expected) const {
    const auto& values = int_descriptor(field);
    require_length(field, values.size());
    return {[values = &values, op, expected](std::size_t atom) {
      if (op == "==") return (*values)[atom] == expected;
      if (op == "!=") return (*values)[atom] != expected;
      if (op == "<") return (*values)[atom] < expected;
      if (op == "<=") return (*values)[atom] <= expected;
      if (op == ">") return (*values)[atom] > expected;
      return (*values)[atom] >= expected;
    }};
  }

  const std::vector<std::string>& string_descriptor(
      const std::string& field) const {
    if (field == "record") return molecule_.record();
    if (field == "name") return molecule_.name();
    if (field == "loc") return molecule_.loc();
    if (field == "resname") return molecule_.resname();
    if (field == "chain") return molecule_.chain();
    if (field == "rescode") return molecule_.rescode();
    if (field == "occupancy") return molecule_.occupancy();
    if (field == "beta") return molecule_.beta();
    if (field == "segname") return molecule_.segname();
    if (field == "element") return molecule_.element();
    if (field == "charge") return molecule_.charge();
    if (field == "charmm_type") return molecule_.charmm_type();
    return molecule_.moltype();
  }

  const std::vector<int>& int_descriptor(const std::string& field) const {
    if (field == "resid") return molecule_.resid();
    if (field == "index") return molecule_.index();
    if (field == "original_index") return molecule_.original_index();
    if (field == "residue_flag") return molecule_.residue_flag();
    return molecule_.original_resid();
  }

  void require_length(const std::string& field, std::size_t actual) const {
    if (actual != molecule_.natoms()) {
      throw std::invalid_argument(selection_error("descriptor '" + field +
                                                  "' length does not match natoms"));
    }
  }

  const Molecule& molecule_;
  std::vector<Token> tokens_;
  std::size_t position_{};
};

SelectionResult collect(const Molecule& molecule,
                        const std::function<bool(std::size_t)>& predicate,
                        const std::string& no_match_message) {
  SelectionResult result;
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    if (predicate(atom)) {
      result.indices.push_back(atom);
    }
  }
  if (result.indices.empty()) {
    result.errors.push_back(no_match_message);
  }
  return result;
}

bool is_comparison_token(const std::string& token) {
  return token == "=" || token == "==" || token == "!=" || token == "<" ||
         token == "<=" || token == ">" || token == ">=";
}

bool is_descriptor_token(const std::string& token) {
  const auto lower = lower_ascii(token);
  return lower == "record" || lower == "name" || lower == "loc" ||
         lower == "resname" || lower == "chain" || lower == "rescode" ||
         lower == "occupancy" || lower == "beta" || lower == "segname" ||
         lower == "element" || lower == "charge" || lower == "moltype" ||
         lower == "charmm_type" || lower == "resid" || lower == "index" ||
         lower == "original_index" || lower == "original_resid" ||
         lower == "residue_flag";
}

bool is_integer_descriptor_token(const std::string& token) {
  const auto lower = lower_ascii(token);
  return lower == "resid" || lower == "index" || lower == "original_index" ||
         lower == "original_resid" || lower == "residue_flag";
}

bool looks_like_python_expression(const std::string& basis) {
  return basis.find("[i]") != std::string::npos;
}

std::string quote_string_literal(const std::string& value) {
  if (value.find('"') == std::string::npos) {
    return "\"" + value + "\"";
  }
  if (value.find('\'') == std::string::npos) {
    return "'" + value + "'";
  }
  throw std::invalid_argument(
      "SASSIE basis string values cannot contain both quote types");
}

std::string comparison_token(std::string token) {
  if (token == "=") {
    return "==";
  }
  return token;
}

std::vector<std::string> tokenize_sassie_basis(const std::string& basis) {
  std::vector<std::string> tokens;
  for (std::size_t offset = 0; offset < basis.size();) {
    const auto value = basis[offset];
    if (std::isspace(static_cast<unsigned char>(value))) {
      ++offset;
      continue;
    }
    if (value == '(' || value == ')') {
      tokens.push_back(std::string(1, value));
      ++offset;
      continue;
    }
    if (value == '"' || value == '\'') {
      const auto quote = value;
      const auto start = ++offset;
      while (offset < basis.size() && basis[offset] != quote) {
        ++offset;
      }
      if (offset == basis.size()) {
        throw std::invalid_argument("unterminated SASSIE basis string literal");
      }
      tokens.push_back(basis.substr(start, offset - start));
      ++offset;
      continue;
    }
    if (value == '=' || value == '!' || value == '<' || value == '>') {
      const auto start = offset++;
      if (offset < basis.size() && basis[offset] == '=') {
        ++offset;
      }
      tokens.push_back(basis.substr(start, offset - start));
      continue;
    }
    const auto start = offset++;
    while (offset < basis.size() &&
           !std::isspace(static_cast<unsigned char>(basis[offset])) &&
           basis[offset] != '(' && basis[offset] != ')' && basis[offset] != '=' &&
           basis[offset] != '!' && basis[offset] != '<' && basis[offset] != '>') {
      ++offset;
    }
    tokens.push_back(basis.substr(start, offset - start));
  }
  return tokens;
}

std::string join_or_name_equals(const std::vector<std::string>& names) {
  std::string expression;
  for (std::size_t index = 0; index < names.size(); ++index) {
    if (index != 0) {
      expression += " or ";
    }
    expression += "name[i] == " + quote_string_literal(names[index]);
  }
  return "(" + expression + ")";
}

BasisExpressionResult alias_expression(const std::string& basis,
                                       SassieBasisContext context) {
  const auto normalized = lower_ascii(trim(basis));
  if (normalized == "all") {
    return {"not name[i] == None", {}};
  }
  if (normalized == "heavy") {
    return {"not name[i][0] == \"H\"", {}};
  }
  if (normalized == "calpha") {
    return {"name[i] == \"CA\"", {}};
  }
  if (normalized == "backbone") {
    if (context == SassieBasisContext::protein) {
      return {join_or_name_equals({"N", "CA", "C"}), {}};
    }
    if (context == SassieBasisContext::nucleic) {
      return {join_or_name_equals({"P", "O5'", "C5'", "C4'", "C3'", "O3'"}),
              {}};
    }
    if (context == SassieBasisContext::nucleic_overlap) {
      return {join_or_name_equals(
                  {"O2P", "O1P", "P", "O5'", "C5'", "C4'", "C3'", "O3'"}),
              {}};
    }
    return {{}, {"SASSIE backbone basis requires an explicit molecule context"}};
  }
  return {{}, {"not an alias"}};
}

std::string value_literal(const std::string& field, const std::string& value) {
  if (value == "None" || value == "True" || value == "False") {
    return value;
  }
  if (is_integer_descriptor_token(field) && parse_int(value)) {
    return value;
  }
  return quote_string_literal(value);
}

std::string translate_sassie_tokens(const std::vector<std::string>& tokens) {
  std::string expression;
  const auto append = [&expression](const std::string& text) {
    if (!expression.empty()) {
      expression += ' ';
    }
    expression += text;
  };
  for (std::size_t index = 0; index < tokens.size();) {
    const auto token = tokens[index];
    const auto lower = lower_ascii(token);
    if (token == "(" || token == ")") {
      append(token);
      ++index;
      continue;
    }
    if (lower == "and" || lower == "or" || lower == "not") {
      append(lower);
      ++index;
      continue;
    }
    if (is_comparison_token(token)) {
      throw std::invalid_argument("unexpected comparison operator '" + token + "'");
    }

    std::string field = "name";
    std::string value = token;
    std::string op = "==";
    if (is_descriptor_token(token)) {
      field = lower;
      if (index + 1 >= tokens.size()) {
        throw std::invalid_argument("descriptor '" + token + "' has no value");
      }
      if (is_comparison_token(tokens[index + 1])) {
        op = comparison_token(tokens[index + 1]);
        if (index + 2 >= tokens.size()) {
          throw std::invalid_argument("descriptor '" + token + "' has no comparison value");
        }
        value = tokens[index + 2];
        index += 3;
      } else {
        value = tokens[index + 1];
        index += 2;
      }
    } else {
      ++index;
    }

    append(field + "[i] " + op + " " + value_literal(field, value));
  }
  return expression;
}

std::vector<std::string> split_commas(const std::string& value) {
  std::vector<std::string> parts;
  std::stringstream stream(value);
  std::string part;
  while (std::getline(stream, part, ',')) {
    parts.push_back(trim(part));
  }
  return parts;
}

}  // namespace

SelectionResult indices_all(const Molecule& molecule) {
  SelectionResult result;
  result.indices.reserve(molecule.natoms());
  for (std::size_t atom = 0; atom < molecule.natoms(); ++atom) {
    result.indices.push_back(atom);
  }
  if (result.indices.empty()) {
    result.errors.push_back("selection found no atoms");
  }
  return result;
}

SelectionResult indices_by_name(const Molecule& molecule, const std::string& name) {
  if (molecule.name().size() != molecule.natoms()) {
    return {{}, {"descriptor 'name' length does not match natoms"}};
  }
  return collect(molecule, [&](std::size_t atom) { return molecule.name()[atom] == name; },
                 "selection found no atoms with name '" + name + "'");
}

SelectionResult indices_by_resname(const Molecule& molecule,
                                   const std::string& resname) {
  if (molecule.resname().size() != molecule.natoms()) {
    return {{}, {"descriptor 'resname' length does not match natoms"}};
  }
  return collect(
      molecule, [&](std::size_t atom) { return molecule.resname()[atom] == resname; },
      "selection found no atoms with resname '" + resname + "'");
}

SelectionResult indices_by_resid_range(const Molecule& molecule, int first_resid,
                                       int last_resid) {
  if (first_resid > last_resid) {
    return {{}, {"resid range lower bound is greater than upper bound"}};
  }
  if (molecule.resid().size() != molecule.natoms()) {
    return {{}, {"descriptor 'resid' length does not match natoms"}};
  }
  return collect(molecule,
                 [&](std::size_t atom) {
                   return molecule.resid()[atom] >= first_resid &&
                          molecule.resid()[atom] <= last_resid;
                 },
                 "selection found no atoms in resid range");
}

std::optional<std::string> basis_expression(const std::string& basis_name) {
  const auto normalized = lower_ascii(basis_name);
  if (normalized == "all") {
    return "not name[i] == None";
  }
  if (normalized == "heavy") {
    return "not name[i][0] == \"H\"";
  }
  return std::nullopt;
}

BasisExpressionResult sassie_basis_expression(const std::string& basis,
                                              SassieBasisContext context) {
  try {
    const auto trimmed = trim(basis);
    if (trimmed.empty()) {
      return {{}, {"SASSIE basis string is empty"}};
    }
    const auto alias = alias_expression(trimmed, context);
    if (alias.ok()) {
      return alias;
    }
    if (lower_ascii(trimmed) == "backbone") {
      return alias;
    }
    if (looks_like_python_expression(trimmed)) {
      return {trimmed, {}};
    }
    return {translate_sassie_tokens(tokenize_sassie_basis(trimmed)), {}};
  } catch (const std::exception& error) {
    return {{}, {"SASSIE basis parse failed: " + std::string(error.what())}};
  }
}

SegmentBasisExpressionResult sassie_segment_basis_expressions(
    const std::string& basis, const std::vector<std::string>& segnames,
    const std::vector<SassieBasisContext>& contexts) {
  if (!contexts.empty() && contexts.size() != segnames.size()) {
    return {{}, {"SASSIE basis contexts must match segment count"}};
  }
  const auto parts = split_commas(basis);
  if (parts.empty() || segnames.empty()) {
    return {{}, {"SASSIE segment basis requires basis and segment names"}};
  }

  const bool single_keyword =
      parts.size() == 1 &&
      (lower_ascii(parts[0]) == "all" || lower_ascii(parts[0]) == "heavy" ||
       lower_ascii(parts[0]) == "backbone");
  if (!single_keyword && parts.size() != segnames.size()) {
    return {{}, {"SASSIE comma basis count must match segment count"}};
  }

  SegmentBasisExpressionResult result;
  result.expressions.reserve(segnames.size());
  for (std::size_t index = 0; index < segnames.size(); ++index) {
    const auto context = contexts.empty() ? SassieBasisContext::generic : contexts[index];
    const auto& basis_part = single_keyword ? parts[0] : parts[index];
    auto expression = sassie_basis_expression(basis_part, context);
    if (!expression.ok()) {
      result.expressions.clear();
      result.errors = expression.errors;
      return result;
    }
    result.expressions.push_back("segname[i] == " +
                                 quote_string_literal(segnames[index]) +
                                 " and ( " + expression.expression + " )");
  }
  return result;
}

SelectionResult select_named_basis(const Molecule& molecule,
                                   const std::string& basis_name) {
  const auto expression = basis_expression(basis_name);
  if (!expression) {
    return {{}, {"unsupported named basis: " + basis_name}};
  }
  return select_indices(molecule, *expression);
}

SelectionResult select_sassie_basis(const Molecule& molecule,
                                    const std::string& basis,
                                    SassieBasisContext context) {
  const auto expression = sassie_basis_expression(basis, context);
  if (!expression.ok()) {
    return {{}, expression.errors};
  }
  return select_indices(molecule, expression.expression);
}

SelectionResult select_indices(const Molecule& molecule,
                               const std::string& expression) {
  try {
    if (expression == "all") {
      return indices_all(molecule);
    }
    Parser parser(molecule, tokenize(expression));
    auto predicate = parser.parse();
    return collect(molecule, predicate.eval,
                   "selection found no atoms using expression '" + expression + "'");
  } catch (const std::exception& error) {
    return {{}, {error.what()}};
  }
}

MaskSelectionResult mask_from_indices(const Molecule& molecule,
                                      const std::vector<std::size_t>& indices) {
  MaskSelectionResult result;
  if (indices.empty()) {
    result.errors.push_back("selection indices are empty");
    return result;
  }
  result.mask.assign(molecule.natoms(), 0);
  for (const auto atom : indices) {
    if (atom >= molecule.natoms()) {
      result.mask.clear();
      result.errors.push_back("selection atom index is out of range");
      return result;
    }
    result.mask[atom] = 1;
  }
  return result;
}

MaskSelectionResult select_mask(const Molecule& molecule,
                                const std::string& expression) {
  const auto selected = select_indices(molecule, expression);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return mask_from_indices(molecule, selected.indices);
}

MaskSelectionResult select_named_basis_mask(const Molecule& molecule,
                                            const std::string& basis_name) {
  const auto selected = select_named_basis(molecule, basis_name);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return mask_from_indices(molecule, selected.indices);
}

MaskSelectionResult select_sassie_basis_mask(const Molecule& molecule,
                                             const std::string& basis,
                                             SassieBasisContext context) {
  const auto selected = select_sassie_basis(molecule, basis, context);
  if (!selected.ok()) {
    return {{}, selected.errors};
  }
  return mask_from_indices(molecule, selected.indices);
}

}  // namespace sasmol
