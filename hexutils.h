#ifndef HEXUTILS_H
#define HEXUITLS_H

#include <string>
#include <sstream>

const char* hexAlphabet = "0123456789abcdef";

bool checkHexString(const std::string& str)
{
	return str.find_first_not_of("0123456789abcdefABCDEF", 0) == std::string::npos;
}

bool checkBinString(const std::string& str)
{
	return str.find_first_not_of("01", 0) == std::string::npos;
}

uint32_t hexStrToUInt32(const std::string& str)
{
	if (!checkHexString(str) && str.size() > 8) return 0;
	std::istringstream sstream(str);
	uint32_t value;
	sstream >> std::hex >> value;
	return value;
}

std::string uint32ToHexString(uint32_t num)
{
	std::string res;
	while (num > 15)
	{
		res += hexAlphabet[num % 16];
		num >>= 4;
	}

	res += hexAlphabet[num];

	if ((res.size() % 8 != 0))
	{
		for (uint8_t i = 0; i < res.size() % 8; i++)
		{
			res += "0";
		}
	}

	return res;
}

uint32_t binStrToUInt32(const std::string& str)
{
	if (!checkBinString(str) && str.size() > 32) return 0;

	uint32_t res = 0;

	for (size_t i = str.size(); i > 0; i--)
	{
		uint32_t temp = (str[i - 1] == '1');
		res += temp << str.size() - i;
	}

	return res;
}

size_t clamp(const size_t& min, const size_t& max, const size_t& value)
{
	return value > min ? (value < max ? value : max) : min;
}

std::string toReverseBinary(uint32_t num)
{
	std::string res;
	while (num != 0)
	{
		res += std::to_string(num % 2);
		num = num / 2;
	}
	if ((res.size() % 32 != 0))
	{
		for (uint8_t i = 0; i < res.size() % 32; i++)
		{
			res += "0";
		}
	}
	return res;
}

#endif