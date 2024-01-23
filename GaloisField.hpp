#ifndef GALOISFIELD_H
#define GALOISFIELD_H

#include <cstdint>
#include <iostream>
#include <bitset>
#include <array>
#include "bigint.h"

// «г≥дно мого вар≥анту


namespace PolyGaloisField
{
	const uint32_t FIELD_DEG = 251;
	const uint32_t NUMBER_SIZE = 32;
	typedef bigint<NUMBER_SIZE> bitvec_t;

	const bitvec_t gen = bitvec_t::fromBinaryString("100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000010011");

	bitvec_t makeGenerator()
	{
		bitvec_t res;
		res.setIthBit(0, 1);
		res.setIthBit(1, 1);
		res.setIthBit(4, 1);
		res.setIthBit(14, 1);
		res.setIthBit(251, 1);

		return res;
	}

	bitvec_t mod(const bitvec_t& a, const bitvec_t& b)
	{
		bitvec_t res(a);
		while (res.bitLength() >= b.bitLength())
		{
			bitvec_t temp(b);
			size_t pos = res.bitLength() - b.bitLength();
			bitvec_t::longShiftBitsToHigh(temp, pos);
			res = res ^ temp;
		}

		return res;
	}

	bitvec_t getZero()
	{
		return bitvec_t::fromConst(0);
	}

	bitvec_t getOne()
	{
		return bitvec_t::fromConst(1);
	}

	bitvec_t add(const bitvec_t& a, const bitvec_t& b)
	{
		return mod(a ^ b, gen);
	}

	bitvec_t mul(const bitvec_t& a, const bitvec_t& b)
	{
		using namespace std;
		bitvec_t res;
		bitvec_t mask = bitvec_t::fromConst(1);

		//cout << "MUL: \n";

		for (int i = 0; i < a.bitLength(); i++)
		{
			bitvec_t t = a & mask;
			t = t * b;

			res = res ^ t;
			bitvec_t::longShiftBitsToHigh(mask, 1);
		}

		return mod(res, gen);
	}


	bitvec_t square(const bitvec_t& a)
	{
		bitvec_t res;

		for (int i = 0, j = 0; i < a.bitLength(); i++, j += 2)
		{
			res.setIthBit(j, a.getIthBit(i));
		}

		return mod(res, gen);
	}

	bitvec_t trace(const bitvec_t& a)
	{
		bitvec_t res;
		bitvec_t t(a);

		for (int i = 0; i < FIELD_DEG; i++)
		{
			res = add(res, t);
			t = square(t);
		}

		return res;
	}
	bitvec_t power(const bitvec_t& a, const bitvec_t& n)
	{
		bitvec_t res = bitvec_t::fromConst(1);

		for (int i = n.bitLength() - 1; i >= 0 ; i--)			
		{
			if (n.getIthBit(i))
			{
				res = mul(res, a);
			}
			if (i != 0)
			{
				res = square(res);
			}
		}

		return res;
	}

	bitvec_t inverse(const bitvec_t& a)
	{
		bitvec_t n = bitvec_t::fromConst(1);
		bitvec_t::longShiftBitsToHigh(n, FIELD_DEG);
		n = n - bitvec_t::fromConst(2);

		return power(a, n);
	}

	bitvec_t makeBitVecFromBitStr(std::string str)
	{
		return bitvec_t::fromBinaryString(str);
	}

	
}


namespace NormalGaloisField
{
	const uint32_t FIELD_DEG = 251;
	using std::bitset;

	typedef std::bitset<FIELD_DEG> bitset_t;

	std::array<std::array<bool, FIELD_DEG>, FIELD_DEG> multMatr;

	bitset_t add(const bitset_t& a, const bitset_t& b)
	{
		return a ^ b;
	}

	bitset_t circularShiftToRight(bitset_t a)
	{
		return (a >> 1) | (a << FIELD_DEG - 1);
	}
	
	bitset_t circularShiftToLeft(bitset_t a)
	{
		return (a << 1) | (a >> FIELD_DEG - 1);
	}

	bitset_t fromBitString(std::string str)
	{
		bitset_t res;

		std::reverse(str.begin(), str.end());

		for (int i = 0; i < str.size() && i < FIELD_DEG; i++)
		{
			switch (str[i])
			{
			case '0':
				res.reset(i);
				break;
			case '1':
				res.set(i);
				break;
			default:
				std::cout << "ERROR! Incorrect bin string!" << std::endl;
				return res;
			}
		}

		return res;
	}

	int32_t modP(int i)
	{
		uint32_t res = 1;
		uint32_t x = 2;
		uint32_t p = 2 * FIELD_DEG + 1;

		while (i > 0)
		{
			if ((i & 1) == 1)
			{
				res = res * x;
				res = res % p;
			}
			x *= x;
			x = x % p;
			i >>= 1;
		}

		return res;
	}

	int32_t mod(int32_t a, int32_t n)
	{
		if (n < 0) n = -n;

		while (a < 0) a += n;

		return a % n;
	} 

	bool multMatrCellChar(int i, int j)
	{
		const int32_t p = (2 * FIELD_DEG + 1);

		const int32_t a = modP(i);
		const int32_t b = modP(j);

		const bool t1 = mod(a + b, p)  == 1;
		const bool t2 = mod(a - b, p) == 1;
		const bool t3 = mod(-a + b, p) == 1;
		const bool t4 = mod(-a - b, p) == 1;

		return t1 || t2 || t3 || t4;
	}

	void genereteMultMatr()
	{
		for (int i = 0; i < FIELD_DEG; i++)
		{
			for (int j = 0; j < FIELD_DEG; j++)
			{
				multMatr[i][j] = multMatrCellChar(i, j);
			}
		}
	}


	bitset_t multVecOnMatr(bitset_t vec)
	{
		bitset_t res;
		for (int i = 0; i < FIELD_DEG; i++)
		{
			for (int j = 0; j < FIELD_DEG; j++)
			{
				res[i] = res[i] ^ (vec[j] & multMatr[FIELD_DEG - 1 - j][i]);
			}
		}

		return res;
	}
	
	bitset_t mult(bitset_t vec1, bitset_t vec2)
	{
		bitset_t res;

		for (int i = FIELD_DEG - 1; i >= 0; i--)
		{
			bitset_t ab = multVecOnMatr(vec1);

			for (int j = 0; j < FIELD_DEG; j++)
			{
				res[i] = res[i] ^ (ab[j] & vec2[FIELD_DEG - 1 - j]);
			}

			vec1 = circularShiftToLeft(vec1);
			vec2 = circularShiftToLeft(vec2);
		}

		return res;
	}


	std::string matrixToString()
	{
		std::string res;
		for (int i = 0; i < FIELD_DEG; i++)
		{
			for (int j = 0; j < FIELD_DEG; j++)
			{
				res += multMatr[i][j] ? '1' : '0';
				res += " | ";
			}
			res += "\n";
		}

		return res;
	}

	uint32_t trace(bitset_t a)
	{
		uint32_t res = 0;
		for (int i = 0; i < FIELD_DEG; i++)
		{
			res ^= a[i];
		}

		return res;
	}

	bitset_t square(bitset_t a)
	{
		return circularShiftToRight(a);
	}

	bitset_t power(bitset_t a, bitset_t n)
	{
		bitset_t res;
		res.flip();

		for (int i = 0; i < n.size(); i++)
		{
			if (n[i])
			{
				res = mult(res, a);
			}
			a = square(a);
			std::cout << i << std::endl;
		}
		return res;
	}

	bitset_t inverse(const bitset_t& a)
	{
		bitset_t b(a);
		uint32_t k = 1;
		uint32_t mask = 1 << 8;

		while (mask != 0)
		{

			bitset_t c(b);

			for (int j = 0; j < k; j++)
			{
				c = square(c);
			}

			b = mult(b, c);
			k = 2 * k;

			if ((FIELD_DEG & mask) != 0)
			{
				std::cout << "ENTER! MASK: " << mask << std::endl;
				b = mult(square(b), a);
				k++;
			}
			std::cout << (mask) << std::endl;

			mask >>= 1;

		}

		return square(b);
	}

	bitset_t slowinverse(const bitset_t& a)
	{
		bitset_t n;
		n.flip();
		n.reset(FIELD_DEG - 1);
		n.reset(0);

		return power(a, n);
	}
}

#endif