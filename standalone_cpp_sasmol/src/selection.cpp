#include "sasmol/selection.hpp"

#include <cctype>
#include <charconv>
#include <functional>
#include <optional>
#include <stdexcept>
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
    auto left = parse_primary();
    while (match(TokenKind::logical_and)) {
      auto right = parse_primary();
      auto left_eval = std::move(left.eval);
      auto right_eval = std::move(right.eval);
      left.eval = [left_eval, right_eval](std::size_t atom) {
        return left_eval(atom) && right_eval(atom);
      };
    }
    return left;
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
    const auto op = expect(TokenKind::comparison, "comparison operator").text;
    const auto value = current();
    if (value.kind != TokenKind::string_literal && value.kind != TokenKind::integer) {
      throw std::invalid_argument(
          selection_error("comparison value must be a string or integer literal"));
    }
    ++position_;

    if (field == "name" || field == "resname" || field == "chain" ||
        field == "segname" || field == "element" || field == "moltype") {
      if (value.kind != TokenKind::string_literal) {
        throw std::invalid_argument(
            selection_error("descriptor '" + field + "' requires a string literal"));
      }
      return string_comparison(field, op, value.text);
    }
    if (field == "resid" || field == "index" || field == "original_index" ||
        field == "original_resid") {
      if (value.kind != TokenKind::integer) {
        throw std::invalid_argument(
            selection_error("descriptor '" + field + "' requires an integer literal"));
      }
      const auto parsed = parse_int(value.text);
      if (!parsed) {
        throw std::invalid_argument(selection_error("invalid integer literal"));
      }
      return int_comparison(field, op, *parsed);
    }
    throw std::invalid_argument(selection_error("unsupported descriptor '" + field + "'"));
  }

  Predicate string_comparison(const std::string& field, const std::string& op,
                              const std::string& expected) const {
    if (op != "==" && op != "!=") {
      throw std::invalid_argument(
          selection_error("string descriptor '" + field +
                          "' only supports == and !="));
    }
    const auto& values = string_descriptor(field);
    require_length(field, values.size());
    return {[values = &values, op, expected](std::size_t atom) {
      const auto equal = (*values)[atom] == expected;
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
    if (field == "name") return molecule_.name();
    if (field == "resname") return molecule_.resname();
    if (field == "chain") return molecule_.chain();
    if (field == "segname") return molecule_.segname();
    if (field == "element") return molecule_.element();
    return molecule_.moltype();
  }

  const std::vector<int>& int_descriptor(const std::string& field) const {
    if (field == "resid") return molecule_.resid();
    if (field == "index") return molecule_.index();
    if (field == "original_index") return molecule_.original_index();
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

}  // namespace sasmol
