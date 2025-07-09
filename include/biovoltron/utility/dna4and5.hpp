#pragma once

#include <algorithm>
#include <istream>
#include <ostream>
#include "istring.hpp"

namespace biovoltron {
    using ichar = std::int8_t;
    using istring = std::basic_string<ichar>;
    using istring_view = std::basic_string_view<ichar>;
    
    
    class dna4 :public biovoltron::Codec {
    public:
        static auto get_seq(const dna4& s) {
            return s.seq;
        }
        dna4(const istring& str) : seq(str) {

        }
        // for original rev_comp and hash function
        operator istring_view() const {
            return istring_view(seq.data(), seq.size());
        }
        // for other char to 0 dna5 is same as Codec
        constexpr static auto ints_dna4 = [] {
            auto ints = std::array<ichar, 128>{};
            ints.fill(0);
            ints['a'] = ints['A'] = 0;
            ints['c'] = ints['C'] = 1;
            ints['g'] = ints['G'] = 2;
            ints['t'] = ints['T'] = 3;
            return ints;
        }();
        constexpr static auto
        to_int(char c) noexcept {
            return ints_dna4[c];
        }
        
        constexpr static auto chars_dna4 = std::array{'A', 'C', 'G', 'T'};
        constexpr static auto
        to_char(int i) noexcept {
            return chars_dna4[i];
        }
        static auto
        to_string(istring_view seq) {
            auto res = std::string{};
            std::ranges::transform(seq, std::back_inserter(res), to_char);
            return res;
        }
        static auto
        to_istring(std::string_view seq) {
            auto res = istring{};
            std::ranges::transform(seq, std::back_inserter(res), to_int);
            return res;
        }
        // for operation get index and change value 
        struct proxy {
            ichar& ref;
            proxy& operator=(char c) {
                ref = to_int(c);
                return *this;
            }
            operator char() const {
                return to_char(ref);
            }
        };

        proxy operator[](std::size_t index) {
            return proxy{seq[index]};
        }

        char operator[](std::size_t index) const {
            return to_char(seq[index]);
        }
        
    private:
        istring seq;
    };
    class dna5 : public biovoltron::Codec{
    public:
        static auto get_seq(const dna5& s) {
            return s.seq;
        }
        dna5(const istring& str) : seq(str) {

        }
        // for original rev_comp and hash function
        operator istring_view() const {
            return istring_view(seq.data(), seq.size());
        }
        // for operation get index and change value 
        struct proxy {
            ichar& ref;
            proxy& operator=(char c) {
                ref = to_int(c);
                return *this;
            }
            operator char() const {
                return to_char(ref);
            }
        };

        proxy operator[](std::size_t index) {
            return proxy{seq[index]};
        }

        char operator[](std::size_t index) const {
            return to_char(seq[index]);
        }
    private:
        istring seq;
    };
    // for 0123_dna4
    inline auto operator""_dna4(char const * s)
    {
        auto is = istring{};
        for (const auto c : std::string_view{s}) is.push_back(c - '0');
        return is;
    }
    // for "ATCG"_dna4
    inline auto operator""_dna4(char const *s, std::size_t n)
    {
        auto is = istring{};
        for (const auto c : std::string_view{s}) is.push_back(dna4::to_int(c));
        return is;
    }
    // for "ATCG"_dna5
    inline auto operator""_dna5(char const *s, std::size_t n)
    {
        auto is = istring{};
        for (const auto c : std::string_view{s}) is.push_back(dna5::to_int(c));
        return is;
    }
    // for 0123_dna5
    inline auto operator""_dna5(char const * s)
    {
        auto is = istring{};
        for (const auto c : std::string_view{s}) is.push_back(c - '0');
        return is;
    }
} // namespace biovoltron

namespace std {
    // output for debug
    inline auto&
    operator<<(ostream& out, biovoltron::dna5 s) {
    return out <<biovoltron::Codec::to_string(biovoltron::dna5::get_seq(s));
    }
    // output for debug
    inline auto&
    operator<<(ostream& out, biovoltron::dna4 s) {
    return out <<biovoltron::dna4::to_string(biovoltron::dna4::get_seq(s));
    }
} // namespace std