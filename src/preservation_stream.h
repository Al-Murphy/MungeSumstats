#ifndef MSS_PRESERVATION_STREAM_H
#define MSS_PRESERVATION_STREAM_H

#include <Rcpp.h>
#include <algorithm>
#include <charconv>
#include <climits>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

// One decoder is used by both collection modes. All fields contribute to the
// bounds and syntax checks, even when the file API only retains selected ones.
class MssDelimitedParser {
    enum class Mode { fresh, rows, columns };
    enum class State { start, quoted, after_quote, plain };
    Mode mode_ = Mode::fresh;
    State state_ = State::start;
    char delimiter_;
    std::size_t limit_, encoded_ = 0, decoded_ = 0, fields_ = 0, selected_ = 0;
    bool cr_ = false, active_ = false, ready_ = false, finished_ = false, eof_seen_ = false;
    bool all_ = true, capture_ = true;
    int width_ = 0, header_width_ = 0, max_rows_ = 0, byte_budget_ = 0;
    int batch_rows_ = 0;
    double emitted_rows_ = 0;
    std::size_t batch_size_ = 0;
    std::vector<int> keep_;
    std::vector<std::vector<std::string>> columns_;
    std::vector<std::string> row_;
    std::string field_;

    void select_mode(Mode mode) {
        if (mode_ != Mode::fresh && mode_ != mode)
            Rcpp::stop("Cannot mix row and column parser collectors");
        mode_ = mode;
    }

    bool selected_field() const {
        return all_ || (selected_ < keep_.size() &&
                        static_cast<std::size_t>(keep_[selected_]) == fields_ + 1);
    }

    void clear_row() {
        row_.clear();
        encoded_ = decoded_ = fields_ = selected_ = 0;
        active_ = ready_ = false;
        capture_ = selected_field();
    }

    void end_field() {
        if (capture_) {
            row_.push_back(std::move(field_));
            field_.clear();
            ++selected_;
        }
        ++fields_;
        ++decoded_;  // One byte per field in the established decoded-byte budget.
        state_ = State::start;
        capture_ = selected_field();
    }

    void append(char c) {
        if (capture_) field_.push_back(c);
        ++decoded_;
    }

    bool read_row(const Rcpp::RawVector& input, R_xlen_t& cursor, bool eof) {
        if (ready_) return true;
        const Rbyte* bytes = RAW(input);
        R_xlen_t length = input.size();
        while (cursor < length) {
            char c = static_cast<char>(bytes[cursor++]);
            if (cr_) {
                cr_ = false;
                if (c == '\n') continue;
            }
            if (c == '\0') Rcpp::stop("NUL byte in tabular input");
            if (++encoded_ > limit_) Rcpp::stop("Input record exceeds max_record_bytes");
            active_ = true;
            if (state_ == State::quoted) {
                if (c == '"') state_ = State::after_quote;
                else append(c);
            } else if (state_ == State::after_quote && c == '"') {
                append(c);
                state_ = State::quoted;
            } else if (c == delimiter_) {
                end_field();
            } else if (c == '\r' || c == '\n') {
                end_field();
                cr_ = c == '\r';
                return ready_ = true;
            } else if (state_ == State::after_quote) {
                Rcpp::stop("Malformed quote: characters after closing quote");
            } else if (c == '"') {
                if (state_ != State::start) Rcpp::stop("Malformed quote in unquoted field");
                state_ = State::quoted;
            } else {
                append(c);
                state_ = State::plain;
            }
        }
        if (eof) {
            if (state_ == State::quoted) Rcpp::stop("Unterminated quoted field");
            if (active_) {
                end_field();
                return ready_ = true;
            }
        }
        return false;
    }

    static Rcpp::List result(const char* state, R_xlen_t cursor, int rows,
                             const Rcpp::List& values) {
        return Rcpp::List::create(Rcpp::_["state"] = state,
                                  Rcpp::_["offset"] = static_cast<int>(cursor),
                                  Rcpp::_["nrow"] = rows, Rcpp::_["columns"] = values);
    }

    Rcpp::List emit_batch(R_xlen_t cursor) {
        Rcpp::List values(columns_.size());
        for (std::size_t j = 0; j < columns_.size(); ++j) {
            values[j] = Rcpp::wrap(columns_[j]);
            columns_[j].clear();
        }
        int rows = batch_rows_;
        emitted_rows_ += rows;
        batch_rows_ = 0;
        batch_size_ = 0;
        return result("batch", cursor, rows, values);
    }

    void configure(int width, const Rcpp::IntegerVector& keep, int max_rows, int budget) {
        if (width == 0) {
            if (header_width_ || keep.size() || max_rows != 1)
                Rcpp::stop("Invalid header parser configuration");
            return;
        }
        if (header_width_ != width) Rcpp::stop("Body width differs from header");
        if (width_) {
            if (width_ != width || max_rows_ != max_rows || byte_budget_ != budget ||
                static_cast<R_xlen_t>(keep_.size()) != keep.size() ||
                !std::equal(keep_.begin(), keep_.end(), keep.begin()))
                Rcpp::stop("Cannot change column parser configuration");
            return;
        }
        int last = 0;
        for (int index : keep) {
            if (index <= last || index > width)
                Rcpp::stop("Selected source column indexes must be increasing and in range");
            last = index;
        }
        keep_.assign(keep.begin(), keep.end());
        columns_.resize(keep_.size());
        width_ = width;
        max_rows_ = max_rows;
        byte_budget_ = budget;
        all_ = false;
        capture_ = selected_field();
    }

public:
    MssDelimitedParser(char delimiter, std::size_t limit): delimiter_(delimiter), limit_(limit) {}

    Rcpp::List feed(const Rcpp::RawVector& input, bool eof) {
        select_mode(Mode::rows);
        if (finished_) Rcpp::stop("Parser already finished");
        R_xlen_t cursor = 0;
        std::vector<std::vector<std::string>> rows;
        while (read_row(input, cursor, eof)) {
            rows.push_back(std::move(row_));
            clear_row();
        }
        finished_ = eof;
        Rcpp::List output(rows.size());
        for (std::size_t j = 0; j < rows.size(); ++j) output[j] = Rcpp::wrap(rows[j]);
        return output;
    }

    Rcpp::List batch(const Rcpp::RawVector& input, int offset, bool eof, int max_rows,
                     int budget, int width, const Rcpp::IntegerVector& keep) {
        select_mode(Mode::columns);
        if (offset < 0 || offset > input.size() || input.size() > INT_MAX ||
            width < 0 || max_rows < 1 || budget < 1)
            Rcpp::stop("Invalid column parser cursor or limits");
        if ((eof && input.size()) || (eof_seen_ && !eof))
            Rcpp::stop("Column parser EOF requires an empty final input");
        eof_seen_ = eof;
        R_xlen_t cursor = offset;
        if (finished_) return result("done", cursor, 0, Rcpp::List(0));
        configure(width, keep, max_rows, budget);
        while (read_row(input, cursor, eof)) {
            if (!width) {
                if (fields_ > INT_MAX) Rcpp::stop("Too many input header fields");
                Rcpp::List header = Rcpp::List::create(Rcpp::wrap(row_));
                header_width_ = static_cast<int>(fields_);
                clear_row();
                return result("batch", cursor, 1, header);
            }
            if (batch_rows_ && decoded_ > static_cast<std::size_t>(budget) - batch_size_)
                return emit_batch(cursor); // This complete row starts the next batch.
            if (fields_ != static_cast<std::size_t>(width))
                Rcpp::stop("Input row has an unexpected number of fields near source row %.0f", emitted_rows_ + 1);
            for (std::size_t j = 0; j < columns_.size(); ++j)
                columns_[j].push_back(std::move(row_[j]));
            ++batch_rows_;
            batch_size_ += decoded_;
            clear_row();
            if (batch_rows_ == max_rows || batch_size_ >= static_cast<std::size_t>(budget))
                return emit_batch(cursor);
        }
        if (!eof) return result("need_input", cursor, 0, Rcpp::List(0));
        if (batch_rows_) return emit_batch(cursor);
        finished_ = true;
        return result("done", cursor, 0, Rcpp::List(0));
    }
};

// Feed native-encoded blocks into the same R text connection used by
// write.table. This keeps filename handling and text-mode newlines with R.
class MssNormalizedWriter {
    Rcpp::Function write_lines_ = Rcpp::Environment::base_env()["writeLines"];
    SEXP connection_;
    std::size_t capacity_;
    std::string buffer_;

    void flush() {
        if (buffer_.empty()) return;
        Rcpp::Shield<SEXP> text(Rf_ScalarString(Rf_mkCharLenCE(
            buffer_.data(), static_cast<int>(buffer_.size()), CE_BYTES)));
        write_lines_(static_cast<SEXP>(text), connection_, Rcpp::Named("sep") = "", Rcpp::Named("useBytes") = true);
        buffer_.clear();
    }

    void bytes(const char* data, std::size_t count) {
        while (count) {
            std::size_t n = std::min(count, capacity_ - buffer_.size());
            buffer_.append(data, n);
            data += n;
            count -= n;
            if (buffer_.size() == capacity_) flush();
        }
    }

    void character(char c) {
        buffer_.push_back(c);
        if (buffer_.size() == capacity_) flush();
    }

    void quoted(const char* data, std::size_t size) {
        character('"');
        while (size) {
            const char* quote = static_cast<const char*>(std::memchr(data, '"', size));
            if (!quote) {
                bytes(data, size);
                break;
            }
            std::size_t n = quote - data + 1;
            bytes(data, n);
            character('"');
            data += n;
            size -= n;
        }
        character('"');
    }

    void string_value(SEXP value) {
        if (value == NA_STRING) {
            bytes("NA", 2);
            return;
        }
        // An ALTREP accessor may return a CHARSXP that the column does not
        // retain. Keep the bytes live across each allocating writeLines call.
        Rcpp::Shield<SEXP> protected_value(value);
        quoted(CHAR(protected_value), LENGTH(protected_value));
    }

    void number(double value) {
        if (ISNAN(value)) {
            bytes("NA", 2);
        } else if (!std::isfinite(value)) {
            if (value < 0) bytes("\"-Inf\"", 6);
            else bytes("\"Inf\"", 5);
        } else {
            char text[64];
            std::size_t n = 0;
            // %.17g spells integers below 1e17 without an exponent. Preserve
            // the sign of zero, which conversion to an integer would discard.
            if (std::abs(value) < 1e17 && std::trunc(value) == value &&
                    !(value == 0 && std::signbit(value))) {
                auto converted = std::to_chars(text, text + sizeof(text), static_cast<long long>(value));
                n = static_cast<std::size_t>(converted.ptr - text);
            } else {
                int length = std::snprintf(text, sizeof(text), "%.17g", value);
                if (length < 0 || static_cast<std::size_t>(length) >= sizeof(text))
                    Rcpp::stop("Cannot format normalized numeric value");
                n = static_cast<std::size_t>(length);
            }
            character('"');
            bytes(text, n);
            character('"');
        }
    }

    // Re-encoding text connections may need whole multibyte characters per
    // write. Delegate those connections, older locales and Windows runtimes to
    // the established platform formatter until there is platform byte evidence.
    static bool native_locale() {
#ifdef _WIN32
        return false;
#else
        if (std::strcmp(std::localeconv()->decimal_point, ".") != 0) return false;
        Rcpp::Function locale_info = Rcpp::Environment::base_env()["l10n_info"];
        Rcpp::List locale = locale_info();
        return Rcpp::as<bool>(locale["UTF-8"]);
#endif
    }

    static bool needs_translation(SEXP column) {
        R_xlen_t size = XLENGTH(column);
        for (R_xlen_t i = 0; i < size; ++i) {
            SEXP element = STRING_ELT(column, i);
            if (element == NA_STRING) continue;
            cetype_t encoding = Rf_getCharCE(element);
            if (encoding != CE_NATIVE && encoding != CE_UTF8) return true;
        }
        return false;
    }

    static bool can_write_native(const Rcpp::DataFrame& table, bool header, bool native_encoding) {
        if (!native_encoding || !native_locale()) return false;
        // Decide before the first output byte. For marked fields that need
        // conversion, the baseline owns translation of the entire table.
        if (header && needs_translation(table.names())) return false;
        for (R_xlen_t j = 0; j < table.size(); ++j) {
            SEXP value = table[j];
            if (TYPEOF(value) == STRSXP && needs_translation(value)) return false;
        }
        return true;
    }

    void base_write(const Rcpp::DataFrame& table, bool header) {
        Rcpp::Environment base = Rcpp::Environment::base_env();
        Rcpp::Function sprintf_r = base["sprintf"], as_character = base["as.character"];
        Rcpp::Function write_table = Rcpp::Environment::namespace_env("utils")["write.table"];
        Rcpp::List encoded(table.size());
        for (R_xlen_t j = 0; j < table.size(); ++j) {
            SEXP value = table[j];
            if (TYPEOF(value) == REALSXP || TYPEOF(value) == INTSXP) {
                Rcpp::CharacterVector text = sprintf_r("%.17g", value);
                for (R_xlen_t i = 0; i < text.size(); ++i) {
                    bool missing = TYPEOF(value) == REALSXP ? ISNAN(REAL_ELT(value, i)) :
                                                             INTEGER_ELT(value, i) == NA_INTEGER;
                    if (missing) text[i] = NA_STRING;
                }
                encoded[j] = text;
            } else {
                encoded[j] = as_character(value);
            }
        }
        encoded.attr("class") = "data.frame";
        encoded.attr("names") = table.names();
        encoded.attr("row.names") = table.attr("row.names");
        write_table(encoded, connection_, Rcpp::Named("sep") = "\t", Rcpp::Named("quote") = true,
                    Rcpp::Named("row.names") = false, Rcpp::Named("col.names") = header,
                    Rcpp::Named("na") = "NA", Rcpp::Named("qmethod") = "double");
    }

public:
    MssNormalizedWriter(SEXP connection, int block_size): connection_(connection) {
        if (block_size < 1) Rcpp::stop("Normalized writer block size must be positive");
        capacity_ = static_cast<std::size_t>(block_size);
        buffer_.reserve(std::min<std::size_t>(capacity_, 262144));
    }

    void write(const Rcpp::DataFrame& table, bool header, bool native_encoding) {
        if (!can_write_native(table, header, native_encoding)) {
            base_write(table, header);
            return;
        }
        if (header) {
            Rcpp::CharacterVector names = table.names();
            for (R_xlen_t j = 0; j < names.size(); ++j) {
                if (j) character('\t');
                string_value(STRING_ELT(names, j));
            }
            character('\n');
        }
        struct Column { SEXP value; int type; };
        std::vector<Column> columns;
        R_xlen_t rows = table.nrows();
        for (R_xlen_t j = 0; j < table.size(); ++j) {
            SEXP value = table[j];
            int type = TYPEOF(value);
            if (XLENGTH(value) != rows || (type != STRSXP && type != LGLSXP &&
                                          type != INTSXP && type != REALSXP))
                Rcpp::stop("Invalid normalized column type or length");
            columns.push_back({value, type});
        }
        for (R_xlen_t i = 0; i < rows; ++i) {
            if (i % 4096 == 0) Rcpp::checkUserInterrupt();
            for (std::size_t j = 0; j < columns.size(); ++j) {
                if (j) character('\t');
                const auto& col = columns[j];
                if (col.type == STRSXP) {
                    string_value(STRING_ELT(col.value, i));
                } else if (col.type == LGLSXP) {
                    int v = LOGICAL_ELT(col.value, i);
                    if (v == NA_LOGICAL) bytes("NA", 2);
                    else if (v) bytes("\"TRUE\"", 6);
                    else bytes("\"FALSE\"", 7);
                } else {
                    double v;
                    if (col.type == REALSXP) v = REAL_ELT(col.value, i);
                    else {
                        int integer = INTEGER_ELT(col.value, i);
                        v = integer == NA_INTEGER ? NA_REAL : static_cast<double>(integer);
                    }
                    number(v);
                }
            }
            character('\n');
        }
        flush();
    }
};

#endif
