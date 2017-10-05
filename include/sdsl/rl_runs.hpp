#ifndef INCLUDED_SDSL_RL_RUNS
#define INCLUDED_SDSL_RL_RUNS

#include "sdsl_concepts.hpp"
#include "int_vector.hpp"
#include "rank_support_v.hpp"
#include "select_support_mcl.hpp"
#include "wt_helper.hpp"
#include "util.hpp"
#include <algorithm> // for std::swap
#include <stdexcept>
#include <vector>
#include <queue>
#include <utility>

//! Namespace for the succinct data structure library.
namespace sdsl
{

template<class t_wt              = wt_gmr<>,              // Wavelet tree type
	 class t_bitvector   = bit_vector,
	 class t_sparse_bitvector   = sd_vector<>,
         class t_rank        = typename t_bitvector::rank_1_type,
         class t_select      = typename t_bitvector::select_1_type,
         class t_sparse_rank        = typename t_sparse_bitvector::rank_1_type,
         class t_sparse_select      = typename t_sparse_bitvector::select_1_type>
class rl_runs
{
    public:

        typedef int_vector<>::size_type              size_type;
        typedef int_vector<>::value_type             value_type;
        typedef random_access_const_iterator<rl_runs> const_iterator;
        typedef const_iterator                       iterator;
        typedef t_bitvector                          bit_vector_type;
        typedef t_sparse_bitvector                   sparse_bit_vector_type;
        typedef t_rank                               rank_1_type;
        typedef t_select                             select_1_type;
        typedef t_sparse_rank                        sparse_rank_1_type;
        typedef t_sparse_select                      sparse_select_1_type;
        typedef wt_tag                               index_category;
        typedef int_alphabet_tag                     alphabet_category;
        typedef t_wt                                 wavelet_tree_type;

protected:

        /* TODO: Change the data type of m_runs_sym and m_X. Try enc_vector */

        size_type         m_size  = 0;
        size_type         m_sigma = 0;     //<- \f$ |\Sigma| \f$
        size_type         m_runs = 0;
        bit_vector_type   m_mapping;       // mapping symbols to consecutive
                                           // alphabet 
        rank_1_type       m_mapping_rank;  // rank support for the
                                           // mapping bit vector 
        t_wt              m_R_wt;          // Default, Golynski's data structure
        int_vector<>      m_runs_sym;      // Number of runs per symbol.
        int_vector<>      m_X;             // Total number of symbols per run,
  
        sparse_bit_vector_type   m_B;         
        sparse_rank_1_type       m_B_rank;
        sparse_select_1_type     m_B_select1;

        void copy(const rl_runs& rl) {
            m_size          = rl.m_size;
            m_sigma         = rl.m_sigma;
            m_runs          = rl.m_runs;
            m_mapping       = rl.m_mapping;
            m_mapping_rank  = rl.m_mapping_rank;
	    m_R_wt          = rl.m_R_wt;
	    m_runs_sym      = rl.m_runs_sym;
	    m_X             = rl.m_X;
	    m_mapping_rank.set_vector(&m_mapping);

	    m_B             = rl.m_B;
            m_B_rank        = rl.m_B_rank;
            m_B_select1     = rl.m_B_select1;
	    m_B_rank.set_vector(&m_B);
	    m_B_select1.set_vector(&m_B);
        }

    public:

        const size_type&       sigma = m_sigma;         //!< Effective alphabet size of the wavelet tree.

        //! Default constructor
        rl_runs() {
        };

        //! Semi-external constructor
        /*! \param buf         File buffer of the int_vector for which the wt_int should be build.
         *  \param size        Size of the prefix of v, which should be indexed.
ermined automatically.
         *    \par Time complexity
         *
         *    \par Space complexity
         *
         */
        template<uint8_t int_width>
        rl_runs(int_vector_buffer<int_width>& buf, size_type size) :
	  m_size(size) {
	  if (0 == m_size)
	    return;
	  size_type n = buf.size();  // set n

	  if (n < m_size) {
	    throw std::logic_error("n="+util::to_string(n)+" < "+util::to_string(m_size)+"=m_size");
	    return;
	  }
	  m_sigma = 0;
	  
	  value_type x = buf[0];  // variable for the biggest value in buf
	  for (size_type i=1; i < m_size; ++i) {
	    if (buf[i] > x )
		x = buf[i];
	    if(buf[i] != buf[i-1])
	      m_runs++;
	  }
	  x++;
	  m_runs++;
	  bit_vector_type mapping(x,0);
	  int_vector<> R(m_runs,0,buf.width());
	  bit_vector_type B_local(buf.size(),0);
	  int_vector<> runs_size(m_runs,0,32);

	  mapping[buf[0]] = 1;
	  size_type l = 0;
	  for (size_type i=1, j=0; i < m_size; ++i) {
	    mapping[buf[i]] = 1;
	    
	    if(buf[i] != buf[i-1]) {
	      R[j] = buf[i-1];
	      B_local[i-1] = 1;
	      runs_size[j] = i-l;
	      l=i;
	      j++;
	    }
	  }
	  R[m_runs-1] = buf[buf.size()-1];
	  B_local[buf.size()-1] = 1;
	  runs_size[m_runs-1] = buf.size()-l;

	  m_mapping.swap(mapping);
	  util::init_support(m_mapping_rank, &m_mapping);
	  m_sigma = m_mapping_rank(x);

	  int_vector<> runs_sym(sigma,0,32);
	  int_vector<> X(m_runs,0,32);
	  for (size_type i=0; i < R.size(); ++i)
	    runs_sym[m_mapping_rank(R[i])]++;	

	  l = 0;
	  size_type tmp = 0;
	  for (size_type i=0; i < runs_sym.size(); ++i) {
	    tmp = runs_sym[i];
	    runs_sym[i] = l;
	    l += tmp;
	  }

	  int_vector<> runs_idx(runs_sym);
	  
	  for (size_type i=0; i < R.size(); ++i) {
	    size_type sym = m_mapping_rank(R[i]);
	    X[runs_idx[sym]] = runs_size[i];
	    runs_idx[sym]++;
	  }
	  runs_idx.resize(0);
       
	  for (size_type i=0; i < sigma; ++i) {
	    size_type ll = runs_sym[i];
	    size_type ul;
	    
	    (i == sigma-1) ? ul = m_runs-1 : ul = runs_sym[i+1]-1;
	    
	    for (size_type i=ll+1; i <= ul; ++i)
	      X[i] += X[i-1];
	  }

	  wavelet_tree_type tmp_wt;
	  construct_im(tmp_wt, R, 0);

	  sparse_bit_vector_type B(B_local);
	  util::init_support(m_B_rank, &B);
	  util::init_support(m_B_select1, &B);

	  m_B.swap(B);
	  m_X.swap(X);
	  m_runs_sym.swap(runs_sym);
	  m_R_wt.swap(tmp_wt);
        }

        //! Copy constructor
        rl_runs(const rl_runs& rl) {
            copy(rl);
        }

        //! Copy constructor
        rl_runs(rl_runs&& rl) {
            *this = std::move(rl);
        }

        //! Assignment operator
        rl_runs& operator=(const rl_runs rl) {
            if (this != &rl) {
                copy(rl);
            }
            return *this;
        }

        //! Assignment move operator
        rl_runs& operator=(rl_runs&& rl) {
            if (this != &rl) {
                m_size          = rl.m_size;
                m_sigma         = rl.m_sigma;
                m_runs          = rl.m_runs;
		m_mapping       = std::move(rl.m_mapping);
		m_mapping_rank  = std::move(rl.m_mapping_rank);
		m_R_wt          = std::move(rl.m_R_wt);
		m_runs_sym      = std::move(rl.m_runs_sym);
		m_X             = std::move(rl.m_X);
		m_mapping_rank.set_vector(&m_mapping);

		m_B             = std::move(rl.m_B);
		m_B_rank        = std::move(rl.m_B_rank);
		m_B_select1     = std::move(rl.m_B_select1);
		m_B_rank.set_vector(&m_B);
		m_B_select1.set_vector(&m_B);
            }
            return *this;
        }

        //! Swap operator
        void swap(rl_runs& rl) {
            if (this != &rl) {
                std::swap(m_size, rl.m_size);
                std::swap(m_sigma,  rl.m_sigma);
                std::swap(m_runs, rl.m_runs);
                m_mapping.swap(rl.m_mapping);
                util::swap_support(m_mapping_rank, rl.m_mapping_rank,
				   &m_mapping, &(rl.m_mapping));
                m_R_wt.swap(rl.m_R_wt);
                m_runs_sym.swap(rl.m_runs_sym);
                m_X.swap(rl.m_X);

		m_B.swap(rl.m_B);
                util::swap_support(m_B_rank, rl.m_B_rank, &m_B, &(rl.m_B));
                util::swap_support(m_B_select1, rl.m_B_select1, &m_B, &(rl.m_B));
            }
        }

        //! Returns the size of the original vector.
        size_type size()const {
            return m_size;
        }

        //! Returns whether the wavelet tree contains no data.
        bool empty()const {
            return m_size == 0;
        }

        //! Returns the number of runs in the sequence
        size_type runs()const {
            return m_runs;
        }

        //! Recovers the i-th symbol of the original vector.
        /*! \param i The index of the symbol in the original vector.
         *  \returns The i-th symbol of the original vector.
         *  \par Precondition
         *       \f$ i < size() \f$
         */
        value_type operator[](size_type i)const {
            assert(i < size());
	    value_type c = m_R_wt[m_B_rank(i)];
	    return c;
        };

        //! Calculates how many symbols c are in the prefix [0..i-1] of the supported vector.
        /*!
         *  \param i The exclusive index of the prefix range [0..i-1], so \f$i\in[0..size()]\f$.
         *  \param c The symbol to count the occurrences in the prefix.
         *    \returns The number of occurrences of symbol c in the prefix [0..i-1] of the supported vector.
         *  \par Time complexity
         *  
         *  \par Precondition
         *  
         */
        size_type rank(size_type i, value_type c)const {
            assert(i <= size());
	    size_type out = 0;
	    size_type idx = m_B_rank(i);
	    size_type r = m_R_wt.rank(idx,c);

	    if(r > 0) { // There are previous runs with symbol c
	      size_type map_sym = m_mapping_rank(c);
	      size_type offset = m_runs_sym[map_sym];
	      out += m_X[offset + r -1];
	    }
	    // Index i is inside a run of symbols c
	    if(m_R_wt[idx] == c) {
	      if(idx)
		idx = m_B_select1(idx)+1;
	      out += (i-idx);
	    }

	    return out;
        };


        //! Calculates the i-th occurrence of the symbol c in the supported vector.
        /*!
         *  \param i The i-th occurrence.
         *  \param c The symbol c.
         *  \par Time complexity
         *  
         *  \par Precondition
         *       \f$ 1 \leq i \leq rank(size(), c) \f$
         */
        size_type select(size_type i, value_type c)const {
            assert(1 <= i and i <= rank(size(), c));
	    
	    size_type map_sym = m_mapping_rank(c);
	    int32_t ll = m_runs_sym[map_sym];
	    int32_t ul;
	    
	    (map_sym == m_sigma-1) ? ul = m_runs-1 : ul =
	      m_runs_sym[map_sym+1]-1;
	    // Find the closest greater/equal value to i in X_c [Binary search]
	    int32_t m;
	    while(ll <= ul) {
	      m = (ll+ul)/2;
	      if(m_X[m] == i)
		break;
	      else if(m_X[m] < i)
		ll = m+1;
	      else
		ul = m-1;
	    }
	    if(m_X[m] < i)
	      m++;

	    size_type local_m = m - m_runs_sym[map_sym];
	    size_type s = m_R_wt.select(local_m+1,c);

	    size_type b;
	    (s == 0) ? b = 0 : b = m_B_select1(s);
	
	    ll = m_runs_sym[map_sym];
	    (m-1 < ll) ? b+=i : b += i - m_X[m-1];

	    // Special case for the first run
	    if (s == 0) b--;
	    
	    return b;
        };

        //! Returns a const_iterator to the first element.
        const_iterator begin()const {
            return const_iterator(this, 0);
        }

        //! Returns a const_iterator to the element after the last element.
        const_iterator end()const {
            return const_iterator(this, size());
        }


        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr, std::string name="")const {
	  structure_tree_node* child = structure_tree::add_child(v, name,
								 util::class_name(*this));

	  size_type written_bytes = 0;
	  written_bytes += write_member(m_size, out, child, "size");
	  written_bytes += write_member(m_sigma, out, child, "sigma");
	  written_bytes += write_member(m_runs, out, child, "runs");
	  written_bytes += m_mapping.serialize(out, child, "mapping");
	  written_bytes += m_mapping_rank.serialize(out, child, "mapping_rank");
	  written_bytes += m_R_wt.serialize(out, child, "R_wt");
	  written_bytes += m_runs_sym.serialize(out, child, "runs_sym");
	  written_bytes += m_X.serialize(out, child, "X");
	  written_bytes += m_B.serialize(out, child, "B");
	  written_bytes += m_B_rank.serialize(out, child, "B_rank");
	  written_bytes += m_B_select1.serialize(out, child, "B_select1");

	  return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_sigma, in);
            read_member(m_runs, in);
	    m_mapping.load(in);
	    m_mapping_rank.load(in, &m_mapping);
	    m_R_wt.load(in);
	    m_runs_sym.load(in);
	    m_X.load(in);

	    m_B.load(in);
	    m_B_rank.load(in, &m_B);
	    m_B_select1.load(in, &m_B);
        }

        void info() {
	  std::cout << "\tsize: " << m_size << std::endl;
	  std::cout << "\tsigma: " << m_sigma << std::endl;
	  std::cout << "\t# of runs: " << m_runs << std::endl;
	  std::cout << "\tTotal space (B): " << size_in_bytes(*this) << std::endl;
	}
};

}// end namespace sdsl
#endif
