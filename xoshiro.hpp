#ifndef XOSHIRO_HPP
#define XOSHIRO_HPP

#include <fstream>
#include <random>


/**
   \file xoshiro.hpp

   \brief contains implementations of a std::random-style class
   of a xoshiro pRNG.

   Based on the code of David Blackman and Sebastiano Vigna at
   http://xoshiro.di.unimi.it/xoshiro256starstar.c
*/
template <typename result_type, int state_size, result_type default_seed>
class xoshiro_engine {
public:
	explicit xoshiro_engine( result_type val = default_seed )
		: s_{0}
	{
		seed(val);
	}

	explicit xoshiro_engine( result_type state[state_size] )
		: s_{0}
	{
		for( int i = 0; i < state_size; ++i ){
			s_[i] = state[i];
		}
	}


	explicit xoshiro_engine( std::seed_seq &q )
		: s_{0}
	{
		seed(q);
	}

	virtual ~xoshiro_engine(){}


	void seed(result_type val = default_seed)
	{
		std::seed_seq seq{val};
		seed(seq);
	}

	void seed( std::seed_seq &q )
	{
		std::vector<uint32_t> seeds(8);
		q.generate(seeds.begin(), seeds.end());

		for( int i = 0; i < 8; i += 2 ){
			uint64_t part1 = seeds[i];
			uint64_t part2 = seeds[i+1];
			uint64_t combined = (part1 << 32) | part2;
			s_[i/2] = combined;
		}
	}

	virtual result_type min()
	{ return my_min(); }

	virtual result_type max()
	{ return my_max(); }

	result_type operator()()
	{ return next(); }

	void discard( unsigned long long step )
	{
		// jump      jumps 2^128 calls ahead.
		// long_jump jumps 2^192 calls ahead.
		for(unsigned long long i = 0; i < step; ++i ){
			next();
		}
	}



	result_type rotl(const result_type x, int k)
	{
		return (x << k) | (x >> (64-k));
	}

private:
	virtual result_type next() = 0;
	virtual result_type my_min() = 0;
	virtual result_type my_max() = 0;

	void generic_jump( result_type jump_arr[state_size] )
	{
		result_type s_tmp[state_size] = {0,0,0,0};
		for(int i = 0; i < state_size; ++i) {
			for(int b = 0; b < 64; ++b) {
				if (jump_arr[i] & UINT64_C(1) << b) {
					// Unroll this loop plz:
					for(int j = 0; j < state_size; ++j ) {
						s_tmp[j] ^= s_[j];
					}
				}
				next();
			}
		}
		for( int i = 0; i < state_size; ++i ){
			s_[i] = s_tmp[i];
		}
	}

	void jump()
	{
		const uint64_t JUMP[] = { 0x180ec6d33cfd0aba,
		                          0xd5a61266f0c9392c,
		                          0xa9582618e03fc9aa,
		                          0x39abdc4529b1661c };
		generic_jump(JUMP);
	}

	void long_jump()
	{
		const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf,
		                               0xc5004e441c522fb3,
		                               0x77710069854ee241,
		                               0x39109bb02acbe635 };
		generic_jump(LONG_JUMP);

	}

protected:
	result_type s_[state_size];
};



/**
   \brief 256-bit state Xoshiro generator.
*/
class xoshiro256 : public xoshiro_engine<uint_fast64_t, 4, 1337>
{
public:
	typedef uint64_t result_type;

	explicit xoshiro256( uint_fast64_t val = 1337 )
		: xoshiro_engine<uint_fast64_t, 4, 1337>(val)
	{ }

	template <class Sseq>
	explicit xoshiro256( Sseq &q )
		: xoshiro_engine<uint_fast64_t, 4, 1337>(q)
	{}

	virtual ~xoshiro256(){}

	friend std::ostream &operator<<( std::ostream &, const xoshiro256 & );
	friend std::istream &operator>>( std::istream &, xoshiro256 & );

	int write_to_binary( std::ostream &out ) const
	{
		const void *state = s_;
		int write_size = 4*sizeof(uint_fast64_t);
		out.write(static_cast<const char*>(state), write_size);
		if( out ) return write_size;
		else      return 0;
	}

	int read_from_binary( std::istream &in )
	{
		void *state = s_;
		int read_size = 4*sizeof(uint_fast64_t);
		in.read(static_cast<char *>(state), read_size);
		if( in ) return read_size;
		else     return 0;
	}



private:
	virtual uint_fast64_t next() override
	{
		const uint_fast64_t result_starstar = rotl(s_[1]*5, 7) * 9;
		const uint_fast64_t t = s_[1] << 17;

		s_[2] ^= s_[0];
		s_[3] ^= s_[1];
		s_[1] ^= s_[2];
		s_[0] ^= s_[3];

		s_[2] ^= t;

		s_[3] = rotl(s_[3], 45);

		return result_starstar;
	}

	virtual uint_fast64_t my_min() override
	{
		return 0;
	}

	virtual uint_fast64_t my_max() override
	{
		const uint_fast64_t largest = 18446744073709551615UL;
		return largest;
	}
};




// Input/output operators:
inline std::ostream &operator<<( std::ostream &out, const xoshiro256 &rng )
{

	for( std::size_t i = 0; i < 4; ++i ){
		out << rng.s_[i] << ' ';
	}
	return out;
}


inline std::istream &operator>>( std::istream &in, xoshiro256 &rng )
{
	for( std::size_t i = 0; i < 4; ++i ){
		in >> rng.s_[i];
	}
	return in;
}



#endif // XOSHIRO_HPP
