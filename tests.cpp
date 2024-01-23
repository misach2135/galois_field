#include <iostream>

#include "GaloisField.hpp"

void testPolyGaloisField()
{
	using namespace std;
	using namespace PolyGaloisField;

	bitvec_t a = bitvec_t::fromHexString("004F23A4027366AFFFEE7219750C7F1FFBC9872D05C4FF1284C879EDC38F868A");
	bitvec_t b = bitvec_t::fromHexString("010125D48680070F1E1930908511C798DC206CE131E6DAD8F7FA4BFA6E76692C");
	bitvec_t n = bitvec_t::fromHexString("077BB0E95DFB7FCAE0EA564F89B3BD431D54FC18D305353B257F839F01EBE608");

	cout << "A: " << a.toBinaryString() << endl;
	cout << "B: " << b.toBinaryString() << endl;

	cout << "A + B: " << add(a, b).toHexString() << endl;
	cout << "A * B: " << (mul(a, b)).toHexString() << endl;
	cout << "A * A: " << (square(a)).toHexString() << endl;
	cout << "A^B: " << (power(a, n)).toHexString() << endl;
	cout << "Tr(A): " << (trace(a)).toHexString() << endl;
	//cout << "A * B: " << (a.multiplyMod(b, gen)).toBinaryString() << endl;
	cout << "A^-1: " << (inverse(a)).toHexString() << endl;
}

void testNormalGaloisField()
{
	using namespace std;
	using namespace NormalGaloisField;

	genereteMultMatr();

	bitset_t t1 = fromBitString("10000111011000110011101001010110011110001000100100000010010110111010100011001100100110011110110100110111011111010101010110000111000000101010011011011110011110010100011010101011100010010010101100101000010000011001010010111010010010000001110110011001101");
	bitset_t t2 = fromBitString("10101110010111011100001011010100111010010111100101010011001110101101001110010111100010000010001100110010010110110011111101110010001000110010010101000011111000011100000000011000100111000110110010111111110001001111001000110101001100001101111101001001011");
	bitset_t n = fromBitString("11100010010100000000011110000000000000110111110011001001111011010010001010011010101010010010111110110101001011111001001111001001101001001010101101110011001001000100011101001010100101110011000001000101110000100011111010010010011000011010000010100001010");

	genereteMultMatr();

	cout << "Add: " << add(t1, t2) << endl;
	cout << "Mult: " << mult(t1, t2) << endl;
	cout << "Trace: " << trace(t1) << endl;
	cout << "Square: " << square(t1) << endl;
	//cout << "Power: "  << power(t1, n) << endl;
	cout << "Inverse: " << inverse(t1) << endl;
}

int main() 
{
	//testNormalGaloisField();
	testPolyGaloisField();

	return 0;
}