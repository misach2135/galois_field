#ifndef BIGINT_H
#define BIGINT_H

#include <iostream>
#include <cstdint>
#include <string>
#include <tuple>
#include <algorithm>

#include "hexutils.h"

template<size_t length>
class bigint
{
template<size_t _len>
friend class bigint;

private:
	static const uint8_t bits = 32;
	uint32_t* digits;

public:

	bigint()
	{
		digits = new uint32_t[length]{0};
		for (size_t i = 0; i < length; i++)
		{
			digits[i] = 0;
		}
	}

	bigint(const bigint<length>& obj) : bigint()
	{
		for (size_t i = 0; i < length; i++)
		{
			digits[i] = obj.digits[i];
		}
	}
	
	~bigint()
	{
		digits = nullptr;
		delete[] digits;
	}

	static bigint<length> fromConst(const uint64_t c)
	{
		bigint<length> res;

		res.digits[0] = (c & UINT32_MAX);
		res.digits[1] = (c >> 32);

		return res;
	}

	static bigint<length> fromHexString(const std::string& hexString)
	{
		if (!checkHexString(hexString)) return bigint();

		bigint<length> res;
		const size_t hexStringSize = hexString.size();
		int64_t i = hexStringSize;

		while(i > 0)
		{
			const size_t digitIndex = (hexStringSize - i) / 8;
			const uint32_t offset = clamp(0, 8, i);
			res.digits[digitIndex] = hexStrToUInt32(hexString.substr(i - offset, offset));
			i -= offset;
		}

		return res;

	}

	static bigint<length> fromBinaryString(const std::string& bin)
	{
		if (!checkBinString(bin)) return bigint();

		bigint<length> res;
		const size_t binStringSize = bin.size();
		int64_t i = binStringSize;

		while (i > 0)
		{
			const size_t digitIndex = (binStringSize - i) / 32;
			const uint32_t offset = clamp(0, 32, i);
			res.digits[digitIndex] = binStrToUInt32(bin.substr(i - offset, offset));
			i -= offset;
		}

		return res;
	}

	static bigint<length> fromHalfLength(const bigint<length / 2>& halfLengthBigInt)
	{
		bigint<length> res;

		for (size_t i = 0; i < length / 2; i++)
		{
			res.digits[i] = halfLengthBigInt.digits[i];
		}

		return res;
	}

	static bigint<length> toHalfLength(const bigint<2 * length>& doubleLengthBigInt)
	{
		bigint<length> res;

		for (size_t i = 0; i < length; i++)
		{
			res.digits[i] = doubleLengthBigInt.digits[i];
		}

		return res;
	}

	static void longShiftBitsToHigh(bigint<length>& number, uint32_t n)
	{
		if (n >= 32)
		{
			longShiftDigitsToHigh(number, n / 32);
			n = n % 32;
		}
		if (n == 0) return;

		uint32_t t = 0;
		const uint32_t mask = (UINT32_MAX) << (32 - n);
		for (size_t i = length - 1; i > 0; --i)
		{
			t = (number.digits[i - 1] & mask) >> (32 - n);
			number.digits[i] = (number.digits[i] << n) | t;
		}
		number.digits[0] <<= n;
	}

	static void longShiftBitsToDown(bigint<length>& number, uint32_t n)
	{
		if (n == 0) return;

		if (n >= 32)
		{
			longShiftDigitsToDown(number, n / 32);
			n = n % 32;
		}

		const uint32_t mask = (UINT32_MAX) >> (32 - n);
		for (size_t i = 0; i < length; i++)
		{
			uint32_t temp = mask & number.digits[i];
			if (i != 0) number.digits[i - 1] = number.digits[i - 1] | (temp << 31);
			number.digits[i] >>= n;
		}
	}

	
	static void longShiftDigitsToHigh(bigint<length>& number, uint32_t n)
	{
		if (n == 0) return;

		for (size_t i = length; i > 0; i--)
		{
			if (i - 1 < n)
			{
				number.digits[i - 1] = 0;
			}
			else
			{
				number.digits[i - 1] = number.digits[i - 1 - n];
			}
		}
	}
	
	static void longShiftDigitsToDown(bigint<length>& number, uint32_t n)
	{
		if (n == 0)
		{
			return;
		}

		for (size_t i = 0; i < length; i++)
		{
			if (i + n >= length)
			{
				number.digits[i] = 0;
				continue;
			}
			number.digits[i] = number.digits[i + n];
		}

	}

	bigint<length> power(const bigint<length>& n) const
	{
		if (n == fromConst(1))
		{
			return *this;
		}

		if (n.isZero())
		{
			return bigint<length>();
		}

		bigint<length> res = fromConst(1);

		for (size_t i = n.bitLength(); i > 0; i--)
		{
			if (n.getIthBit(i - 1))
			{
				res = res * *this;
			}
			if (i - 1 != 0)
			{
				res = res * res;
			}
		}
		return res;
	}

	bool isEven() const
	{
		return (digits[0] & (uint32_t)1) == ((uint32_t)0);
	}

	size_t bitLength() const
	{
		size_t i = length - 1;
		size_t res = 0;
		while (digits[i] == 0)
		{
			i--;
		}
		
		uint32_t n = digits[i];

		while (n != 0)
		{
			n >>= 1;
			res++;
		}
		res += i * 32;
		return res;
	}

	size_t countTrailingZeros() const
	{
		size_t i = 0;
		size_t res = 0;
		while (digits[i] == 0)
		{
			res += 32;
			i++;
		}

		uint32_t digit = digits[i];
		const uint32_t mask = 1;

		while (!(digit & mask))
		{
			digit >>= 1;
			res++;
		}

		return res;
	}

	size_t getNumberLength() const
	{
		size_t i = length - 1;
		size_t res = 0;
		while (digits[i] == 0)
		{
			i--;
		}

		return i + 1;
	}

	bool getIthBit(size_t i) const
	{
		const uint32_t mask = 1 << (i % 32);
		return (digits[i / 32] >> (i % 32)) & (uint32_t)1;
	}

	void setIthBit(size_t i, bool bit)
	{
		const uint32_t mask = 1 << (i % 32);
		switch (bit)
		{
		case 0:
			digits[i / 32] = digits[i / 32] & ~mask;
			return;
		case 1:
			digits[i / 32] = digits[i / 32] | mask;
			return;
		}
	}

	std::pair<bigint<length>, bigint<length> > longDivision(const bigint<length>& a, const bigint<length>& b) const
	{
		const size_t k = b.bitLength();
		bigint<length> r(a);
		bigint<length> q;

		size_t t = 0;

		while (r >= b)
		{
			t = r.bitLength();
			bigint<length> c(b);
			longShiftBitsToHigh(c, t - k);
			if (r < c)
			{
				t--;
				longShiftBitsToDown(c, 1);
			}
			r = r - c;
			q.setIthBit(t - k, 1);
		}

		return std::make_pair(q, r);
	}

	bigint<length> operator+ (const bigint<length>& b) const
	{
		uint32_t carry = 0;
		bigint<length> res;
		for (size_t i = 0; i < length; i++)
		{
			uint64_t temp = (uint64_t)digits[i] + (uint64_t)b.digits[i] + (uint64_t)carry;
			res.digits[i] = temp & UINT32_MAX;
			carry = uint8_t(temp >> bits);
		}

		return res;
	}

	bigint<length> operator- (const bigint<length>& b) const
	{
		bigint<length> res;
		uint8_t borrow = 0;
		for (size_t i = 0; i < length; i++)
		{
			int64_t temp = (uint64_t)digits[i] - (uint64_t)b.digits[i] - (uint64_t)borrow;
			if (temp >= 0)
			{
				res.digits[i] = temp;
				borrow = 0;
			}
			else
			{
				res.digits[i] = temp + ((uint64_t)1 << bits);
				borrow = 1;
			}
		}
		return res;
	}

	bigint<length> operator* (const bigint<length>& b) const
	{
		bigint<length> res;

		for (size_t i = 0; i < length; i++)
		{
			bigint<length> temp = *this * b.digits[i];
			longShiftDigitsToHigh(temp, i);
			res = res + temp;
		}
		return res;
	}
	
	bigint<length> operator* (const uint32_t digit) const
	{
		uint32_t carry = 0;
		bigint<length> res;

		for (size_t i = 0; i < length; i++)
		{
			uint64_t temp = (uint64_t)this->digits[i] * (uint64_t)digit + carry;
			res.digits[i] = temp & UINT32_MAX;
			carry = temp >> bits;
		}

		return res;
	}

	bigint<length> operator/ (const bigint<length>& b) const
	{
		auto res = longDivision(*this, b);
		return res.first;
	}

	bigint<length> operator% (const bigint<length>& n) const
	{
		auto res = longDivision(*this, n);
		return res.second;
	}

	bigint<length> operator& (const bigint<length>& b) const
	{
		bigint<length> res(*this);
		for (size_t i = 0; i < length; i++)
		{
			res.digits[i] &= b.digits[i];
		}
		return res;
	}

	bigint<length> operator| (const bigint<length>& b) const
	{
		bigint<length> res(*this);
		for (size_t i = 0; i < length; i++)
		{
			res.digits[i] |= b.digits[i];
		}
		return res;
	}

	bigint<length> operator^ (const bigint<length>& b) const
	{
		bigint<length> res(*this);
		for (size_t i = 0; i < length; i++)
		{
			res.digits[i] ^= b.digits[i];
		}
		return res;
	}

	bool isZero() const
	{
		size_t i = length;
		while (i > 0)
		{
			if (digits[i - 1] != 0) return false;
			i--;
		}
		return true;
	}

	bool isOne() const
	{
		for (size_t i = length - 1; i > 0; i--)
		{
			if (digits[i] != 0) return false;
		}

		return digits[0] == (uint32_t)1;
	}

	bool isInt64Const(uint64_t _const) const
	{
		for (size_t i = length - 1; i > 1; i--)
		{
			if (digits[i] != 0) return false;
		}

		uint64_t res = ((uint64_t)digits[0] | ((uint64_t)digits[1]) << 32);
		return res == _const;
	}

	friend bool operator== (const bigint<length>& a, const bigint<length>& b)
	{
		size_t i = length;

		while (i != 0 && a.digits[i - 1] == b.digits[i - 1])
		{
			i--;
		}

		return i == 0;
	}

	friend bool operator== (const bigint<length>& a, const uint64_t _const)
	{
		return a.isInt64Const(_const);
	}

	friend bool operator< (const bigint<length>& a, const bigint<length>& b)
	{
		size_t i = length - 1;

		while (a.digits[i] == b.digits[i])
		{
			if (i == 0) break;
			i--;
		}

		return a.digits[i] < b.digits[i];
	}

	friend bool operator> (const bigint<length>& a, const bigint<length>& b)
	{
		return b < a;
	}

	friend bool operator>= (const bigint<length>& a, const bigint<length>& b)
	{
		return !(a < b);
	}

	friend bool operator<= (const bigint<length>& a, const bigint<length>& b)
	{
		return !(a > b);
	}

	std::string toString() const
	{
		std::string res;
		for (size_t i = length; i > 0; i--)
		{
			res += std::to_string(digits[i - 1]);
			switch (i)
			{
			case 1:
				res += ';';
				break;
			default:
				res += '_';
				break;
			}
		}
		return res;
	}

	std::string toHexString() const 
	{
		if (isZero()) return "0";
		std::string res = "";
		const size_t len = getNumberLength();
		for (size_t i = 0; i < len; i++)
		{
			res += uint32ToHexString(digits[i]);
		}

		std::reverse(res.begin(), res.end());

		return res;
	}

	std::string toBinaryString(bool full = false) const
	{
		if (isZero()) return "0";

		std::string res;

		uint32_t last = length;
		while (!full && last != 0 && digits[last - 1] == 0)
		{
			last--;
		}

		for (size_t i = 0; i < last; i++)
		{
			if (digits[i] == 0)
			{
				res += "00000000000000000000000000000000";
				continue;
			}

			res += toReverseBinary(digits[i]);
		}

		std::reverse(res.begin(), res.end());

		return res;
	}

	static bigint<length> gcd(bigint<length> a, bigint<length> b)
	{
		bigint<length> d = fromConst(1);

		while (a.isEven() && b.isEven())
		{
			longShiftBitsToDown(a, 1);
			longShiftBitsToDown(b, 1);
			longShiftBitsToHigh(d, 1);
		}

		while (a.isEven())
		{
			longShiftBitsToDown(a, 1);
		}

		while (!b.isZero())
		{
			while (b.isEven())
			{
				longShiftBitsToDown(b, 1);
			}
			bigint<length> temp = std::max(a, b) - std::min(a, b);
			a = std::min(a, b);
			b = temp;
		}

		return d * a;
	}

	static bigint<length> lcm(const bigint<length>& a, const bigint<length>& b)
	{
		return (a * b) / gcd(a, b);
	}

	bigint<length> addMod(const bigint<length>& b, const bigint<length>& mod) const
	{
		bigint<length> res = *this + b;
		res = res % mod;
		return res;
	}
	
	bigint<length> subMod(const bigint<length>& b, const bigint<length>& mod) const
	{
		bigint<length> res = *this - b;
		res = res % mod;
		return res;
	}

	bigint<length> multiplyMod(const bigint<length>& b, const bigint<length>& mod) const
	{
		bigint<2 * length> res = bigint<length * 2>::fromHalfLength(*this * b);
		res = res % bigint<length*2>::fromHalfLength(mod);
		return toHalfLength(res);
	}

	bigint<length> barretReduction(const bigint<length>& mod)
	{
		if (*this > mod.power(fromConst(2)))
		{
			return *this % mod;
		}

		size_t k = mod.getNumberLength();

		bigint<length> mu;
		mu.setIthBit(2*k + 1, 1);
		std::cout << "mu : " << mu.toBinaryString() << std::endl;
		mu = mu / mod;

		std::cout << "a : " << toBinaryString() << std::endl;
		std::cout << "mod : " << mod.toBinaryString() << std::endl;

		bigint<length> q(*this);
		q.longShiftDigitsToDown(q, k - 1);
		q = q * mu;
		q.longShiftDigitsToDown(q, k + 1);
		bigint<length> r = *this - q * mod;
		return r;
	}

	bigint<length> powerMod(const bigint<length>& n, const bigint<length>& mod) const
	{
		if (n == fromConst(1))
		{
			return *this;
		}

		if (n.isZero())
		{
			return bigint<length>();
		}

		bigint<length> res = fromConst(1);

		for (size_t i = n.bitLength(); i > 0; i--)
		{
			if (n.getIthBit(i - 1))
			{
				res = res * *this;
			}
			if (i - 1 != 0)
			{
				res = res * res;
			}
			res = res % mod;
		}
		return res;
	}

};

#endif