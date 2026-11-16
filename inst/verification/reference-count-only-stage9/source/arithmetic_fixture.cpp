#define main count_cli_main
#include "count_reference.cpp"
#undef main
int main() {
  U64 value = std::numeric_limits<U64>::max();
  try { add(value, 1); return 1; } catch (const std::runtime_error&) {}
  if (value != std::numeric_limits<U64>::max()) return 2;
  CountStream stream;
  stream.buffer.bytes = std::numeric_limits<U64>::max() - 2;
  try { stream << "abc"; return 3; } catch (const std::exception&) {}
  if (stream.buffer.bytes != std::numeric_limits<U64>::max() - 2) return 4;
  return 0;
}
