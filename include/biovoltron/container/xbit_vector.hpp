#pragma once

#include <cassert>
#include <climits>
#include <limits>
#include <stdexcept>

namespace biovoltron::detail {
/**
 * @brief A reference to an N-bit wide field within a block of memory.
 *
 * Provides read/write access to N-bit values packed in a `Block`.
 *
 * @tparam N Bit-width of each field.
 * @tparam Block Unsigned integral type used as storage block.
 */
template<std::size_t N, std::unsigned_integral Block>
class XbitReference {
 public:
  
  /**
   * Bit mask used to isolate N-bit field.
   */
  constexpr static Block mask =
    (std::numeric_limits<Block>::max() >> (sizeof(Block) * CHAR_BIT - N));

 private:
  Block* seg_;
  const std::size_t shift_;

 public:
  /**
   * Construct a reference to the N-bit field at a given offset.
   *
   * @param seg Pointer to the block containing the field.
   * @param offset The N-bit field index within the block.
   */
  constexpr XbitReference(Block* seg, std::size_t offset) noexcept
  : seg_(seg), shift_(offset * N) { }
  /**
   * Implicit conversion to the referenced N-bit value.
   *
   * @return The value of the N-bit field as a uint8_t.
   */
  constexpr operator std::uint8_t() const noexcept {
    return *seg_ >> shift_ & mask;
  }

  /**
   * Prefix increment the field value.
   *
   * @return Reference to the updated field.
   */
  constexpr XbitReference&
  operator++() noexcept {
    return operator=(*this + 1);
  }

  /**
   * Postfix increment the field value.
   *
   * @return Copy of the value before increment.
   */
  constexpr std::uint8_t
  operator++(int) noexcept {
    std::uint8_t tmp = *this;
    ++*this;
    return tmp;
  }

  /**
   * Prefix decrement the field value.
   *
   * @return Reference to the updated field.
   */
  constexpr XbitReference&
  operator--() noexcept {
    return operator=(*this - 1);
  }

  /**
   * Postfix decrement the field value.
   *
   * @return Copy of the value before decrement.
   */
  constexpr std::uint8_t
  operator--(int) noexcept {
    std::uint8_t tmp = *this;
    --*this;
    return tmp;
  }

  /**
   * Assign a new value to the N-bit field.
   *
   * @param x Value to assign.
   * @return Reference to this field.
   */
  constexpr XbitReference&
  operator=(std::uint8_t x) noexcept {
    *seg_ &= ~(mask << shift_);
    *seg_ |= (x & mask) << shift_;
    return *this;
  }

  /**
   * Assign from another `XbitReference`.
   *
   * @param x Another field reference.
   * @return Reference to this field.
   */
  constexpr XbitReference&
  operator=(const XbitReference& x) noexcept {
    return operator=(static_cast<std::uint8_t>(x));
  }

  /**
   * Dummy const assignment (no-op).
   */
  constexpr void
  operator=(std::uint8_t) const noexcept { }
};
/**
 * Swap two XbitReferences.
 *
 * @tparam N Bit-width of each field.
 * @tparam Block Unsigned integral type used as storage block.
 * @param x First field reference.
 * @param y Second field reference.
 */
template<std::size_t N, std::unsigned_integral Block>
inline void
swap(XbitReference<N, Block> x, XbitReference<N, Block> y) noexcept {
  std::uint8_t t = x;
  x = y;
  y = t;
}
/**
 * Swap an XbitReference with a uint8_t.
 *
 * @tparam N Bit-width of each field.
 * @tparam Block Unsigned integral type used as storage block.
 * @param x XbitReference value.
 * @param y uint8_t value.
 */
template<std::size_t N, std::unsigned_integral Block>
inline void
swap(XbitReference<N, Block> x, std::uint8_t& y) noexcept {
  std::uint8_t t = x;
  x = y;
  y = t;
}
/**
 * Swap a uint8_t with an XbitReference.
 *
 * @tparam N Bit-width of each field.
 * @tparam Block Unsigned integral type used as storage block.
 * @param x uint8_t value.
 * @param y XbitReference value.
 */
template<std::size_t N, std::unsigned_integral Block>
inline void
swap(std::uint8_t& x, XbitReference<N, Block> y) noexcept {
  std::uint8_t t = x;
  x = y;
  y = t;
}
/**
 * @brief Base class for Xbit iterators accessing N-bit fields.
 *
 * @tparam N Number of bits per element.
 * @tparam Block Unsigned integral type used as storage block.
 */
template<std::size_t N, std::unsigned_integral Block>
class XbitIteratorBase {
  public:
   /**
    * Iterator category tag.
    */
    using iterator_category = std::random_access_iterator_tag;
    /**
     * Type of the value pointed to by the iterator.
     */
    using value_type = std::uint8_t;
    /**
     * Type for iterator differences.
     */
    using difference_type = std::ptrdiff_t;
    /**
     * Pointer type for the iterator.
     */
    using pointer = void;
    /**
     * Number of N-bit elements stored in one Block.
     */
  constexpr static std::size_t xbits_per_block = sizeof(Block) * CHAR_BIT / N;

  protected:
    /**
     * Pointer to the current block segment.
     */
    Block* seg_;
    /**
     * Offset within the current block segment.
     */
    std::size_t offset_;

  public:
    /**
     * Construct an iterator pointing to a specific block and offset.
     *
     * @param seg Pointer to the current block.
     * @param offset Index of the N-bit element within the block.
     */
    constexpr XbitIteratorBase(Block* seg, std::size_t offset) noexcept
    : seg_(seg), offset_(offset) { }

    /**
     * Compute the distance between two iterators.
     *
     * @param x Left-hand iterator.
     * @param y Right-hand iterator.
     * @return Difference between the two iterators.
     */
    constexpr friend difference_type
    operator-(const XbitIteratorBase& x, const XbitIteratorBase& y) {
      return (x.seg_ - y.seg_) * xbits_per_block + x.offset_ - y.offset_;
    }

    /**
     * Equality operator for two iterators.
     *
     * @param other Another iterator to compare with.
     * @return true if both iterators point to the same position, false otherwise.
     */
    constexpr bool
    operator==(const XbitIteratorBase& other) const noexcept = default;

    /**
     * Three-way comparison operator for iterators.
     *
     * @param other Another iterator to compare with.
     * @return A value representing the comparison result.
     */
    constexpr auto
    operator<=>(const XbitIteratorBase& other) const noexcept {
      if (auto cmp = seg_ <=> other.seg_; cmp != 0)
        return cmp;
      return offset_ <=> other.offset_;
    }

  protected:
    /**
     * Increment the iterator position.
     *
     * Moves to the next N-bit element.
     */
    constexpr void
    bump_up() {
      if (offset_ != xbits_per_block - 1)
        ++offset_;
      else {
        offset_ = 0;
        ++seg_;
      }
    }

    /**
     * Decrement the iterator position.
     *
     * Moves to the previous N-bit element.
     */
    constexpr void
    bump_down() {
      if (offset_ != 0)
        --offset_;
      else {
        offset_ = xbits_per_block - 1;
        --seg_;
      }
    }

    /**
     * Increment the iterator position by n elements.
     *
     * @param n Number of elements to increment.
     */
    constexpr void
    incr(difference_type n) {
      if (n >= 0)
        seg_ += (n + offset_) / xbits_per_block;
      else
        seg_ += static_cast<difference_type>(n - xbits_per_block + offset_ + 1)
              / static_cast<difference_type>(xbits_per_block);
      n &= (xbits_per_block - 1);
      offset_ = (n + offset_) & (xbits_per_block - 1);
    }
};
/**
 * @brief Random-access iterator for mutable N-bit elements.
 * 
 * @tparam N Number of bits per element.
 * @tparam Block Unsigned integral type used as storage block.
 */
template<std::size_t N, std::unsigned_integral Block>
struct XbitIterator : public XbitIteratorBase<N, Block> {
  using Base = XbitIteratorBase<N, Block>;
  using iterator_category = Base::iterator_category;
  using value_type = Base::value_type;
  using difference_type = Base::difference_type;
  using pointer = Base::pointer;
  using reference = XbitReference<N, Block>;
  using iterator = XbitIterator;

  template<std::size_t, std::unsigned_integral>
  friend class XbitConstIterator;
  /**
   * Default-construct a null iterator.
   */
  constexpr XbitIterator() noexcept : XbitIteratorBase<N, Block>(nullptr, 0) { }
  /**
   * Construct an iterator at a specific block and offset.
   * @param seg Pointer to block containing element.
   * @param offset Index of N-bit element within block.
   */
  constexpr XbitIterator(Block* seg, std::size_t offset) noexcept
  : XbitIteratorBase<N, Block>(seg, offset) { }

  /**
   * Dereference to access the current N-bit element.
   * @return A reference to the N-bit element.
   */
  constexpr reference
  operator*() const noexcept {
    return reference(this->seg_, this->offset_);
  }

  /**
   * Access the N-bit element at a specific index.
   * @param n Element offset from current position.
   * @return A reference to the N-bit element at index n.
   */
  constexpr reference
  operator[](difference_type n) const {
    return *(*this + n);
  }

  /**
   * Pre-increment the iterator.
   * @return Reference to the incremented iterator.
   */
  constexpr iterator&
  operator++() {
    this->bump_up();
    return *this;
  }

  /**
   * Post-increment the iterator.
   * @return A copy of the iterator before incrementing.
   */
  constexpr iterator
  operator++(int) {
    iterator tmp = *this;
    this->bump_up();
    return tmp;
  }

  /**
   * Pre-decrement the iterator.
   * @return Reference to the decremented iterator.
   */
  constexpr iterator&
  operator--() {
    this->bump_down();
    return *this;
  }

  /**
   * Post-decrement the iterator.
   * @return A copy of the iterator before decrementing.
   */
  constexpr iterator
  operator--(int) {
    iterator tmp = *this;
    this->bump_down();
    return tmp;
  }

  /**
   * Increment the iterator position by n elements.
   * @param n Number of elements to increment.
   */
  constexpr iterator&
  operator+=(difference_type n) {
    this->incr(n);
    return *this;
  }

  /**
   * Decrement the iterator position by n elements.
   * @param n Number of elements to decrement.
   */
  constexpr iterator&
  operator-=(difference_type n) {
    return *this += -n;
  }

  /**
   * Add an offset to the iterator.
   * @param n Number of elements to add.
   * @return A new iterator positioned at the sum of the current position and n.
   */
  constexpr iterator
  operator+(difference_type n) const {
    iterator tmp(*this);
    tmp += n;
    return tmp;
  }

  /**
   * Subtract an offset from the iterator.
   * @param n Number of elements to subtract.
   * @return A new iterator positioned at the difference of the current position and n.
   */
  constexpr iterator
  operator-(difference_type n) const {
    iterator tmp(*this);
    tmp -= n;
    return tmp;
  }

  /**
   * Add an offset to the iterator.
   * @param n Number of elements to add.
   * @param it The iterator to advance.
   * @return A new iterator positioned at the sum of the current position and n.
   */
  constexpr friend iterator
  operator+(difference_type n, const iterator& it) {
    return it + n;
  }
};

/**
 * @brief Random-access iterator for const N-bit elements.
 *
 * @tparam N Number of bits per element.
 * @tparam Block Unsigned integral type used as storage block.
 */
template<std::size_t N, std::unsigned_integral Block>
struct XbitConstIterator : public XbitIteratorBase<N, Block> {
  using Base = XbitIteratorBase<N, Block>;
  using iterator_category = Base::iterator_category;
  using value_type = Base::value_type;
  using difference_type = Base::difference_type;
  using pointer = Base::pointer;
  using reference = value_type;
  using const_reference = value_type;
  using const_iterator = XbitConstIterator;

  /**
   * Default-construct a null const iterator.
   */
  constexpr XbitConstIterator() noexcept
  : XbitIteratorBase<N, Block>(nullptr, 0) { }

  /**
   * Construct a const iterator at a specific block and offset.
   * 
   * @param seg Pointer to block containing element.
   * @param offset Index of N-bit element within block.
   */
  constexpr XbitConstIterator(Block* seg, std::size_t offset) noexcept
  : XbitIteratorBase<N, Block>(seg, offset) { }

  /**
   * Construct from a mutable iterator.
   * 
   * @param x Mutable iterator to copy position from.
   */
  
   constexpr XbitConstIterator(const XbitIterator<N, Block>& x) noexcept
  : XbitIteratorBase<N, Block>(x.seg_, x.offset_) { }
  
  /**
   * Dereference to read the current N-bit element.
   * 
   * @return Value of the current element.
   */
  constexpr const_reference
  operator*() const noexcept {
    return XbitReference<N, Block>(this->seg_, this->offset_);
  }

  /**
   * Access element at relative position.
   * 
   * @param n Element offset from current position.
   * @return Value of the element.
   */
  constexpr const_reference
  operator[](difference_type n) const {
    return *(*this + n);
  }

  constexpr const_iterator&
  operator++() {
    this->bump_up();
    return *this;
  }

  constexpr const_iterator
  operator++(int) {
    const_iterator tmp = *this;
    this->bump_up();
    return tmp;
  }

  constexpr const_iterator&
  operator--() {
    this->bump_down();
    return *this;
  }

  constexpr const_iterator
  operator--(int) {
    const_iterator tmp = *this;
    this->bump_down();
    return tmp;
  }

  constexpr const_iterator&
  operator+=(difference_type n) {
    this->incr(n);
    return *this;
  }

  constexpr const_iterator&
  operator-=(difference_type n) {
    return *this += -n;
  }

  constexpr const_iterator
  operator+(difference_type n) const {
    const_iterator tmp(*this);
    tmp += n;
    return tmp;
  }

  constexpr const_iterator
  operator-(difference_type n) const {
    const_iterator tmp(*this);
    tmp -= n;
    return tmp;
  }

  constexpr friend const_iterator
  operator+(difference_type n, const const_iterator& it) {
    return it + n;
  }
};

/**
 * @brief Base class for XbitVector providing common exception helpers.
 *
 * Defines utility methods for throwing standard exceptions in runtime
 * (suppressed during constant evaluation).
 */
class XbitVectorBase {
 protected:
  
  /**
   * Default constructor.
   */
  constexpr XbitVectorBase() = default;
  
  /**
   * Throw std::length_error if not evaluated at compile time.
   *
   * @throw std::length_error Always throws at runtime.
   */
  constexpr void
  throw_length_error() const {
    if (!std::is_constant_evaluated())
      throw std::length_error("XbitVector");
  }

  /**
   * Throw std::out_of_range if not evaluated at compile time.
   *
   * @throw std::out_of_range Always throws at runtime.
   */
  constexpr void
  throw_out_of_range() const {
    if (!std::is_constant_evaluated())
      throw std::out_of_range("XbitVector");
  }
};

/**
 * @ingroup container
 * 
 * @brief Fixed-bit-width packed vector container.
 * 
 * @tparam N Bits per element.
 * @tparam Block Unsigned integral type used as storage block.
 * @tparam Allocator Allocator type for block storage.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
  requires(!std::same_as<Block, bool>)
class XbitVector : private detail::XbitVectorBase {
 public:
  using value_type = std::uint8_t;
  using block_type = Block;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = detail::XbitReference<N, block_type>;
  using const_reference = value_type;
  using iterator = detail::XbitIterator<N, block_type>;
  using const_iterator = detail::XbitConstIterator<N, block_type>;
  using pointer = iterator;
  using const_pointer = const value_type*;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;

  constexpr static std::size_t xbits_per_block = iterator::xbits_per_block;

 private:
  block_type* begin_{};
  size_type size_{};
  size_type cap_{};
  /*[[no_unique_address]]*/ allocator_type alloc_{};

 public:

  /**
   * Default constructor.
   */
  constexpr XbitVector() noexcept(
    std::is_nothrow_default_constructible_v<allocator_type>);

  /**
   * Constructor with allocator.
   */
  constexpr explicit XbitVector(const allocator_type& a) noexcept;

  /**
   * Destructor.
   */
  constexpr ~XbitVector();

  /**
   * Constructors for various initializations.
   */
  constexpr explicit XbitVector(size_type n);

  constexpr explicit XbitVector(size_type n, const allocator_type& a);

  constexpr XbitVector(size_type n, const value_type& x);

  constexpr XbitVector(size_type n, const value_type& x,
                       const allocator_type& a);

  constexpr XbitVector(std::input_iterator auto first,
                       std::input_iterator auto last);

  constexpr XbitVector(std::input_iterator auto first,
                       std::input_iterator auto last, const allocator_type& a);

  constexpr XbitVector(std::forward_iterator auto first,
                       std::forward_iterator auto last);

  constexpr XbitVector(std::forward_iterator auto first,
                       std::forward_iterator auto last,
                       const allocator_type& a);

  constexpr XbitVector(const XbitVector& v);

  constexpr XbitVector(const XbitVector& v, const allocator_type& a);

  /**
   * Copy assignment operator.
   *
   * @param v Source vector to copy from.
   * @return Reference to this vector.
   */
  constexpr XbitVector&
  operator=(const XbitVector& v);

  constexpr XbitVector(std::initializer_list<value_type> il);

  constexpr XbitVector(std::initializer_list<value_type> il,
                       const allocator_type& a);

  /**
   * Move constructor.
   *
   * @param v Source vector to move from.
   */
  constexpr XbitVector(XbitVector&& v) noexcept;

  constexpr XbitVector(XbitVector&& v, const allocator_type& a);

  /**
   * Move assignment operator.
   *
   * @param v Source vector to move from.
   * @return Reference to this vector.
   */
  constexpr XbitVector&
  operator=(XbitVector&& v) noexcept(
    allocator_traits::propagate_on_container_move_assignment::value
    || allocator_traits::is_always_equal::value);

  /**
   * Replaces contents with elements from initializer list.
   *
   * @param il Initializer list to copy from.
   * @return Reference to this vector.
   */
  constexpr XbitVector&
  operator=(std::initializer_list<value_type> il) {
    assign(il.begin(), il.end());
    return *this;
  }

  /**
   * Assigns new contents to the vector.
   */
  constexpr void
  assign(std::input_iterator auto first, std::input_iterator auto last);

  constexpr void
  assign(std::forward_iterator auto first, std::forward_iterator auto last);

  constexpr void
  assign(size_type n, const value_type& x);

  constexpr void
  assign(std::initializer_list<value_type> il) {
    assign(il.begin(), il.end());
  }

  /**
   * Returns the allocator used by the vector.
   *
   * @return Allocator used by the vector.
   */
  constexpr allocator_type
  get_allocator() const noexcept {
    return allocator_type(this->alloc_);
  }

  constexpr size_type
  max_size() const noexcept;

  /**
   * Returns the number of elements that can be held without reallocation.
   */
  constexpr size_type
  capacity() const noexcept {
    return internal_cap_to_external(cap_);
  }

  /**
   * Returns the number of stored elements.
   */
  constexpr size_type
  size() const noexcept {
    return size_;
  }

  /**
   * Returns the number of storage blocks used.
   */
  constexpr size_type
  num_blocks() const noexcept {
    return empty() ? 0 : external_cap_to_internal(size());
  }

  /**
   * Checks if the vector is empty.
   *
   * @return true if the vector contains no elements, false otherwise.
   */
  [[nodiscard]] constexpr bool
  empty() const noexcept {
    return size_ == 0;
  }

  constexpr void
  reserve(size_type n);

  constexpr void
  shrink_to_fit() noexcept;

  /** @name Iterators */
  /// @{
  constexpr iterator
  begin() noexcept {
    return make_iter(0);
  }

  constexpr const_iterator
  begin() const noexcept {
    return make_iter(0);
  }

  constexpr iterator
  end() noexcept {
    return make_iter(size_);
  }

  constexpr const_iterator
  end() const noexcept {
    return make_iter(size_);
  }

  constexpr reverse_iterator
  rbegin() noexcept {
    return reverse_iterator(end());
  }

  constexpr const_reverse_iterator
  rbegin() const noexcept {
    return const_reverse_iterator(end());
  }

  constexpr reverse_iterator
  rend() noexcept {
    return reverse_iterator(begin());
  }

  constexpr const_reverse_iterator
  rend() const noexcept {
    return const_reverse_iterator(begin());
  }

  constexpr const_iterator
  cbegin() const noexcept {
    return make_iter(0);
  }

  constexpr const_iterator
  cend() const noexcept {
    return make_iter(size_);
  }

  constexpr const_reverse_iterator
  crbegin() const noexcept {
    return rbegin();
  }

  constexpr const_reverse_iterator
  crend() const noexcept {
    return rend();
  }
  ///@}

  /** @name Element access */
  ///@{
  constexpr reference
  operator[](size_type n) {
    return *make_iter(n);
  }

  constexpr const_reference
  operator[](size_type n) const {
    return *make_iter(n);
  }

  constexpr reference
  at(size_type n);

  constexpr const_reference
  at(size_type n) const;

  constexpr reference
  front() {
    return *begin();
  }

  constexpr const_reference
  front() const {
    return *begin();
  }

  constexpr reference
  back() {
    return *(end() - 1);
  }

  constexpr const_reference
  back() const {
    return *(end() - 1);
  }

  constexpr block_type*
  data() noexcept {
    return begin_;
  }

  constexpr const block_type*
  data() const noexcept {
    return begin_;
  }
  ///@}

  /** @name Modifiers */
  ///@{
  constexpr void
  push_back(const value_type& x);

  template<typename... Args>
  constexpr reference
  emplace_back(Args&&... args) {
    push_back(value_type(std::forward<Args>(args)...));
    return this->back();
  }

  constexpr void
  pop_back() {
    --size_;
  }

  template<typename... Args>
  constexpr iterator
  emplace(const_iterator position, Args&&... args) {
    return insert(position, value_type(std::forward<Args>(args)...));
  }

  constexpr iterator
  insert(const_iterator position, const value_type& x);

  constexpr iterator
  insert(const_iterator position, size_type n, const value_type& x);

  constexpr iterator
  insert(const_iterator position, std::input_iterator auto first,
         std::input_iterator auto last);

  constexpr iterator
  insert(const_iterator position, std::forward_iterator auto first,
         std::forward_iterator auto last);

  constexpr iterator
  insert(const_iterator position, std::initializer_list<value_type> il) {
    return insert(position, il.begin(), il.end());
  }

  constexpr iterator
  erase(const_iterator position);

  constexpr iterator
  erase(const_iterator first, const_iterator last);

  constexpr void
  clear() noexcept {
    size_ = 0;
  }

  constexpr void
  swap(XbitVector&) noexcept;

  constexpr void
  resize(size_type sz, value_type x = 0);

  constexpr void
  flip() noexcept;
  ///@}

  /** @name Comparisons */
  ///@{
  constexpr bool
  operator==(const XbitVector& other) const {
    return size() == other.size() && std::equal(begin(), end(), other.begin());
  }

  constexpr auto
  operator<=>(const XbitVector& other) const {
    return std::lexicographical_compare_three_way(begin(), end(), other.begin(),
                                                  other.end());
  }
  ///@}

 private:
  constexpr bool
  invariants() const;

  constexpr void
  invalidate_all_iterators();

  constexpr void
  vallocate(size_type n);

  constexpr void
  vdeallocate() noexcept;
/**
 * @brief Converts an internal capacity (in storage blocks) to an external capacity (in elements).
 *
 * @param n Internal capacity, expressed as the number of allocated blocks.
 * @return External capacity, expressed as the number of storable elements.
 * @see external_cap_to_internal()
 */
  constexpr static size_type
  internal_cap_to_external(size_type n) noexcept {
    return n * xbits_per_block;
  }
/**
 * @brief Converts an external capacity (in elements) to an internal capacity (in storage blocks).
 *
 * @param n External capacity, expressed as the number of elements to store.
 * @return Internal capacity, expressed as the number of required blocks.
 *
 * @see internal_cap_to_external()
 */
  constexpr static size_type
  external_cap_to_internal(size_type n) noexcept {
    return (n - 1) / xbits_per_block + 1;
  }

  constexpr static size_type
  align_it(size_type new_size) noexcept {
    return (new_size + (xbits_per_block - 1)) / xbits_per_block
         * xbits_per_block;
  }

  constexpr size_type
  recommend(size_type new_size) const;

  constexpr void
  construct_at_end(size_type n, value_type x);

  constexpr void
  construct_at_end(std::forward_iterator auto first,
                   std::forward_iterator auto last);

  /**
   * @brief Creates a mutable iterator pointing to the specified element position.
   *
   * @param pos Zero-based logical element index.
   * @return Iterator referring to the specified element.
   */
  constexpr iterator
  make_iter(size_type pos) noexcept {
    return iterator(begin_ + pos / xbits_per_block,
                    pos & (xbits_per_block - 1));
  }

  /**
   * @brief Creates a constant iterator pointing to the specified element position.
   * 
   * @param pos Zero-based logical element index.
   * @return Constant iterator referring to the specified element.
   */
  constexpr const_iterator
  make_iter(size_type pos) const noexcept {
    return const_iterator(begin_ + pos / xbits_per_block,
                          pos & (xbits_per_block - 1));
  }

  constexpr iterator
  const_iterator_cast(const_iterator p) noexcept {
    return begin() + (p - cbegin());
  }

  constexpr void
  copy_assign_alloc(const XbitVector& v) {
    if constexpr (allocator_traits::propagate_on_container_copy_assignment::
                    value) {
      if (alloc_ != v.alloc_)
        vdeallocate();
      alloc_ = v.alloc_;
    }
  }

  constexpr void
  move_assign(XbitVector& v) noexcept(
    std::is_nothrow_move_assignable_v<allocator_type>);

  constexpr void
  move_assign_alloc(XbitVector& v) noexcept(
    !allocator_traits::propagate_on_container_move_assignment::value
    || std::is_nothrow_move_assignable_v<allocator_type>) {
    if constexpr (allocator_traits::propagate_on_container_move_assignment::
                    value)
      alloc_ = std::move(v.alloc_);
  }
};

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::invalidate_all_iterators() { }

/**
 * @brief Allocates internal storage for the container.
 *
 * @param n External capacity (number of elements) to reserve.
 *
 * @throw std::length_error If `n` exceeds `max_size()`.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::vallocate(size_type n) {
  if (n > max_size())
    this->throw_length_error();
  n = external_cap_to_internal(n);
  this->begin_ = allocator_traits::allocate(this->alloc_, n);
  this->size_ = 0;
  this->cap_ = n;
}

/**
 * @brief Deallocates all internal storage and resets the container.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::vdeallocate() noexcept {
  if (this->begin_ != nullptr) {
    allocator_traits::deallocate(this->alloc_, this->begin_, this->cap_);
    invalidate_all_iterators();
    this->begin_ = nullptr;
    this->size_ = this->cap_ = 0;
  }
}

/**
 * @brief Returns the maximum number of elements the container can hold.
 *
 * @return Maximum possible element count.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::size_type
XbitVector<N, Block, Allocator>::max_size() const noexcept {
  size_type amax = allocator_traits::max_size(alloc_);
  size_type nmax = std::numeric_limits<size_type>::max() / 2;
  if (nmax / xbits_per_block <= amax)
    return nmax;
  return internal_cap_to_external(amax);
}

/**
 * @brief Suggests a new capacity for growth.
 *
 * @param new_size Desired size after growth.
 * @return Recommended new capacity.
 *
 * @throw std::length_error If `new_size` exceeds `max_size()`.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::size_type
XbitVector<N, Block, Allocator>::recommend(size_type new_size) const {
  const size_type ms = max_size();
  if (new_size > ms)
    this->throw_length_error();
  const size_type cap = capacity();
  if (cap >= ms / 2)
    return ms;
  return std::max(2 * cap, align_it(new_size));
}

/**
 * @brief Constructs `n` elements with the given value at the end.
 *
 * Expands the container size by `n` and fills the newly created
 * positions with the provided value `x`.
 *
 * @param n Number of elements to append.
 * @param x Value to fill in the new elements.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::construct_at_end(size_type n, value_type x) {
  size_type old_size = this->size_;
  this->size_ += n;
  std::fill_n(make_iter(old_size), n, x);
}
/**
 * @brief Constructs elements from a range at the end of the container.
 *
 * Expands the container size by the number of elements in the given
 * iterator range `[first, last)`, and copies the elements into the
 * new positions.
 *
 * @tparam ForwardIt A forward iterator type.
 * @param first Iterator to the first element to copy.
 * @param last  Iterator past the last element to copy.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::construct_at_end(
  std::forward_iterator auto first, std::forward_iterator auto last) {
  size_type old_size = this->size_;
  this->size_ += std::distance(first, last);
  std::copy(first, last, make_iter(old_size));
}

/**
 * @brief Default constructor.
 *
 * Constructs an empty XbitVector with zero capacity and size.
 * The allocator is default-constructed.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector() noexcept(
  std::is_nothrow_default_constructible_v<allocator_type>) { }

/**
 * @brief Constructs an empty XbitVector with a specific allocator.
 *
 * @param a Allocator to use for memory allocation.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  const allocator_type& a) noexcept
: cap_(0), alloc_(a) { }

/**
 * @brief Constructs a XbitVector with n elements, each initialized to zero.
 *
 * @param n Number of elements to create.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(size_type n) {
  if (n > 0) {
    vallocate(n);
    construct_at_end(n, 0);
  }
}

/**
 * @brief Constructs a XbitVector with n elements (value 0) and a specific allocator.
 *
 * @param n Number of elements to create.
 * @param a Allocator to use.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(size_type n,
                                                      const allocator_type& a)
: cap_(0), alloc_(a) {
  if (n > 0) {
    vallocate(n);
    construct_at_end(n, 0);
  }
}
/**
 * @brief Constructs a XbitVector with n elements, each initialized to x.
 *
 * @param n Number of elements.
 * @param x Value to initialize each element with.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(size_type n,
                                                      const value_type& x) {
  if (n > 0) {
    vallocate(n);
    construct_at_end(n, x);
  }
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(size_type n,
                                                      const value_type& x,
                                                      const allocator_type& a)
: cap_(0), alloc_(a) {
  if (n > 0) {
    vallocate(n);
    construct_at_end(n, x);
  }
}

/**
 * @brief Constructs a XbitVector from an input iterator range [first, last).
 *
 * @tparam InputIt Input iterator type.
 * @param first Iterator to the first element to copy.
 * @param last  Iterator past the last element to copy.
 *
 * @throw Any exception thrown by element copy or memory allocation.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::input_iterator auto first, std::input_iterator auto last) {
  try {
    for (; first != last; ++first) push_back(*first);
  } catch (...) {
    if (begin_ != nullptr)
      allocator_traits::deallocate(alloc_, begin_, cap_);
    invalidate_all_iterators();
    throw;
  }
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::input_iterator auto first, std::input_iterator auto last,
  const allocator_type& a)
: cap_(0), alloc_(a) {
  try {
    for (; first != last; ++first) push_back(*first);
  } catch (...) {
    if (begin_ != nullptr)
      allocator_traits::deallocate(alloc_, begin_, cap_);
    invalidate_all_iterators();
    throw;
  }
}

/**
 * @brief Constructs a XbitVector from a forward iterator range [first, last).
 *
 * @tparam ForwardIt Forward iterator type.
 * @param first Iterator to the first element.
 * @param last  Iterator past the last element.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::forward_iterator auto first, std::forward_iterator auto last) {
  const size_type n = std::distance(first, last);
  if (n > 0) {
    vallocate(n);
    construct_at_end(first, last);
  }
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::forward_iterator auto first, std::forward_iterator auto last,
  const allocator_type& a)
: cap_(0), alloc_(a) {
  const size_type n = std::distance(first, last);
  if (n > 0) {
    vallocate(n);
    construct_at_end(first, last);
  }
}

/**
 * @brief Constructs a XbitVector from an initializer list.
 *
 * @param il Initializer list of elements.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::initializer_list<value_type> il) {
  const size_type n = il.size();
  if (n > 0) {
    vallocate(n);
    construct_at_end(il.begin(), il.end());
  }
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(
  std::initializer_list<value_type> il, const allocator_type& a)
: cap_(0), alloc_(a) {
  const size_type n = il.size();
  if (n > 0) {
    vallocate(n);
    construct_at_end(il.begin(), il.end());
  }
}

/**
 * @brief Destructor.
 *
 * Releases all allocated storage and invalidates all iterators.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::~XbitVector() {
  if (begin_ != nullptr)
    allocator_traits::deallocate(alloc_, begin_, cap_);
  invalidate_all_iterators();
}

/**
 * @brief Copy constructor.
 *
 * @param v Source vector to copy.
 *
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(const XbitVector& v)
: cap_(0),
  alloc_(allocator_traits::select_on_container_copy_construction(v.alloc_)) {
  if (v.size() > 0) {
    vallocate(v.size());
    construct_at_end(v.begin(), v.end());
  }
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(const XbitVector& v,
                                                      const allocator_type& a)
: cap_(0), alloc_(a) {
  if (v.size() > 0) {
    vallocate(v.size());
    construct_at_end(v.begin(), v.end());
  }
}

/**
 * @brief Copy assignment operator.
 *
 * Replaces the contents of this vector with a copy of another.
 *
 * @param v Source vector to copy from.
 * @return Reference to this vector.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>&
XbitVector<N, Block, Allocator>::operator=(const XbitVector& v) {
  if (this != &v) {
    copy_assign_alloc(v);
    if (v.size_) {
      if (v.size_ > capacity()) {
        vdeallocate();
        vallocate(v.size_);
      }
      std::copy(v.begin_, v.begin_ + external_cap_to_internal(v.size_), begin_);
    }
    size_ = v.size_;
  }
  return *this;
}

/**
 * @brief Move constructor.
 *
 * Moves the contents of another XbitVector into this one.
 *
 * @param v Source vector to move.
 */

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(XbitVector&& v) noexcept
: begin_(v.begin_), size_(v.size_), cap_(v.cap_), alloc_(v.alloc_) {
  v.begin_ = nullptr;
  v.size_ = 0;
  v.cap_ = 0;
}

/**
 * @brief Move constructor with specified allocator.
 *
 * Moves the contents if the allocator matches; otherwise, copies elements.
 *
 * @param v Source vector to move.
 * @param a Allocator to use.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::XbitVector(XbitVector&& v,
                                                      const allocator_type& a)
: cap_(0), alloc_(a) {
  if (a == allocator_type(v.alloc_)) {
    this->begin_ = v.begin_;
    this->size_ = v.size_;
    this->cap_ = v.cap_;
    v.begin_ = nullptr;
    v.cap_ = v.size_ = 0;
  } else if (v.size() > 0) {
    vallocate(v.size());
    construct_at_end(v.begin(), v.end());
  }
}

/**
 * @brief Move assignment operator.
 *
 * Moves the contents from another vector or copies if allocators differ.
 *
 * @param v Source vector to move.
 * @return Reference to this vector.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>&
XbitVector<N, Block, Allocator>::operator=(XbitVector&& v) noexcept(
  allocator_traits::propagate_on_container_move_assignment::value
  || allocator_traits::is_always_equal::value) {
  if constexpr (allocator_traits::propagate_on_container_move_assignment::value)
    move_assign(v);
  else {
    if (alloc_ != v.alloc_)
      assign(v.begin(), v.end());
    else
      move_assign(v);
  }
  return *this;
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::move_assign(XbitVector& c) noexcept(
  std::is_nothrow_move_assignable_v<allocator_type>) {
  vdeallocate();
  move_assign_alloc(c);
  this->begin_ = c.begin_;
  this->size_ = c.size_;
  this->cap_ = c.cap_;
  c.begin_ = nullptr;
  c.cap_ = c.size_ = 0;
}

/**
 * @brief Assigns n copies of value x to the vector.
 *
 * @param n Number of elements.
 * @param x Value to fill.
 */

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::assign(size_type n, const value_type& x) {
  size_ = 0;
  if (n > 0) {
    size_type c = capacity();
    if (n <= c)
      size_ = n;
    else {
      XbitVector v(alloc_);
      v.reserve(recommend(n));
      v.size_ = n;
      swap(v);
    }
    std::fill_n(begin(), n, x);
  }
  invalidate_all_iterators();
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::assign(std::input_iterator auto first,
                                        std::input_iterator auto last) {
  clear();
  for (; first != last; ++first) push_back(*first);
}

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::assign(std::forward_iterator auto first,
                                        std::forward_iterator auto last) {
  clear();
  difference_type ns = std::distance(first, last);
  assert(ns >= 0 && "invalid range specified");
  const size_type n = ns;
  if (n) {
    if (n > capacity()) {
      vdeallocate();
      vallocate(n);
    }
    construct_at_end(first, last);
  }
}

/**
 * @brief Requests a change in capacity to at least n elements.
 *
 * If n exceeds current capacity, storage is reallocated and elements moved.
 *
 * @param n Minimum capacity requested.
 */

template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::reserve(size_type n) {
  if (n > capacity()) {
    XbitVector v(this->alloc_);
    v.vallocate(n);
    v.construct_at_end(this->begin(), this->end());
    swap(v);
    invalidate_all_iterators();
  }
}

/**
 * @brief Reduces the capacity of the vector to fit its current size.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::shrink_to_fit() noexcept {
  if (external_cap_to_internal(size()) > cap_) {
    try {
      XbitVector(*this, allocator_type(alloc_)).swap(*this);
    } catch (...) { }
  }
}

/**
 * @brief Provides checked access to an element at index `n`.
 * 
 * @param n Index of the element.
 * @return Reference to the element at position `n`.
 * 
 * @throws std::out_of_range if n >= size().
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::reference
XbitVector<N, Block, Allocator>::at(size_type n) {
  if (n >= size())
    this->throw_out_of_range();
  return (*this)[n];
}
/**
 * @brief Provides checked access to an element at index `n` (const version).
 *
 * @param n Index of the element.
 * @return Const reference to the element at position `n`.
 * 
 * @throws std::out_of_range if n >= size().
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::const_reference
XbitVector<N, Block, Allocator>::at(size_type n) const {
  if (n >= size())
    this->throw_out_of_range();
  return (*this)[n];
}

/**
 * @brief Appends an element to the end of the vector.
 *
 * @param x Element to add.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::push_back(const value_type& x) {
  if (this->size_ == this->capacity())
    reserve(recommend(this->size_ + 1));
  ++this->size_;
  back() = x;
}

/**
 * @brief Inserts a single element at the specified position.
 *
 * @param position Iterator pointing to the insertion position.
 * @param x Element to insert.
 * @return Iterator pointing to the inserted element.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::insert(const_iterator position,
                                        const value_type& x) {
  iterator r;
  if (size() < capacity()) {
    const_iterator old_end = end();
    ++size_;
    std::copy_backward(position, old_end, end());
    r = const_iterator_cast(position);
  } else {
    XbitVector v(alloc_);
    v.reserve(recommend(size_ + 1));
    v.size_ = size_ + 1;
    r = std::copy(cbegin(), position, v.begin());
    std::copy_backward(position, cend(), v.end());
    swap(v);
  }
  *r = x;
  return r;
}

/**
 * @brief Inserts `n` copies of `x` at the specified position.
 *
 * @param position Iterator pointing to the insertion position.
 * @param n Number of copies to insert.
 * @param x Element to insert.
 * @return Iterator pointing to the first of the newly inserted elements.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::insert(const_iterator position, size_type n,
                                        const value_type& x) {
  iterator r;
  size_type c = capacity();
  if (n <= c && size() <= c - n) {
    const_iterator old_end = end();
    size_ += n;
    std::copy_backward(position, old_end, end());
    r = const_iterator_cast(position);
  } else {
    XbitVector v(alloc_);
    v.reserve(recommend(size_ + n));
    v.size_ = size_ + n;
    r = std::copy(cbegin(), position, v.begin());
    std::copy_backward(position, cend(), v.end());
    swap(v);
  }
  std::fill_n(r, n, x);
  return r;
}

/**
 * @brief Inserts a range of elements [first, last) at the specified position
 *        (input iterator version).
 *
 * @param position Iterator pointing to the insertion position.
 * @param first Start of the input range.
 * @param last End of the input range.
 * @return Iterator pointing to the first of the newly inserted elements.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr typename XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::insert(const_iterator position,
                                        std::input_iterator auto first,
                                        std::input_iterator auto last) {
  difference_type off = position - begin();
  iterator p = const_iterator_cast(position);
  iterator old_end = end();
  for (; size() != capacity() && first != last; ++first) {
    ++this->size_;
    back() = *first;
  }
  XbitVector v(alloc_);
  if (first != last) {
    try {
      v.assign(first, last);
      difference_type old_size = old_end - begin();
      difference_type old_p = p - begin();
      reserve(recommend(size() + v.size()));
      p = begin() + old_p;
      old_end = begin() + old_size;
    } catch (...) {
      erase(old_end, end());
      throw;
    }
  }
  p = std::rotate(p, old_end, end());
  insert(p, v.begin(), v.end());
  return begin() + off;
}
/**
 * @brief Inserts a range of elements [first, last) at the specified position
 *        (forward iterator version).
 *
 * @param position Iterator pointing to the insertion position.
 * @param first Start of the forward range.
 * @param last End of the forward range.
 * @return Iterator pointing to the first of the newly inserted elements.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::insert(const_iterator position,
                                        std::forward_iterator auto first,
                                        std::forward_iterator auto last) {
  const difference_type n_signed = std::distance(first, last);
  assert(n_signed >= 0 && "invalid range specified");
  const size_type n = n_signed;
  iterator r;
  size_type c = capacity();
  if (n <= c && size() <= c - n) {
    const_iterator old_end = end();
    size_ += n;
    std::copy_backward(position, old_end, end());
    r = const_iterator_cast(position);
  } else {
    XbitVector v(alloc_);
    v.reserve(recommend(size_ + n));
    v.size_ = size_ + n;
    r = std::copy(cbegin(), position, v.begin());
    std::copy_backward(position, cend(), v.end());
    swap(v);
  }
  std::copy(first, last, r);
  return r;
}
/**
 * @brief Removes the element at the specified position.
 *
 * @param position Iterator pointing to the element to remove.
 * @return Iterator pointing to the element following the erased one.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::erase(const_iterator position) {
  iterator r = const_iterator_cast(position);
  std::copy(position + 1, this->cend(), r);
  --size_;
  return r;
}

/**
 * @brief Removes elements in the range [first, last).
 *
 * @param first Start of the range to remove.
 * @param last End of the range to remove.
 * @return Iterator pointing to the element following the last removed one.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr XbitVector<N, Block, Allocator>::iterator
XbitVector<N, Block, Allocator>::erase(const_iterator first,
                                       const_iterator last) {
  iterator r = const_iterator_cast(first);
  difference_type d = last - first;
  std::copy(last, this->cend(), r);
  size_ -= d;
  return r;
}

/**
 * @brief Swaps the contents of this vector with another.
 * 
 * @param x Vector to swap with.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::swap(XbitVector& x) noexcept {
  std::swap(this->begin_, x.begin_);
  std::swap(this->size_, x.size_);
  std::swap(this->cap_, x.cap_);
  if constexpr (allocator_traits::propagate_on_container_swap::value)
    std::swap(this->alloc_, x.alloc_);
}

/**
 * @brief Resizes the vector to contain `sz` elements.
 *
 * @param sz New size of the vector.
 * @param x Value to use for new elements if the vector grows.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::resize(size_type sz, value_type x) {
  size_type cs = size();
  if (cs < sz) {
    iterator r;
    size_type c = capacity();
    size_type n = sz - cs;
    if (n <= c && cs <= c - n) {
      r = end();
      size_ += n;
    } else {
      XbitVector v(alloc_);
      v.reserve(recommend(size_ + n));
      v.size_ = size_ + n;
      r = std::copy(cbegin(), cend(), v.begin());
      swap(v);
    }
    std::fill_n(r, n, x);
  } else
    size_ = sz;
}

/**
 * @brief Flips all bits in the vector.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr void
XbitVector<N, Block, Allocator>::flip() noexcept {
  size_type n = 0;
  for (block_type* p = begin_; n < size_; ++p, n += xbits_per_block) *p = ~*p;
}

/**
 * @brief Checks internal container invariants.
 *
 * This function verifies that the internal state of the container
 * is consistent and does not violate basic assumptions:
 * - If `begin_` is `nullptr`, both `size_` and `cap_` must be zero.
 * - If `begin_` is not `nullptr`, `cap_` must be non-zero.
 * - `size_` must never exceed `capacity()`.
 *
 * @return true if internal state is valid, false otherwise.
 */
template<std::size_t N, std::unsigned_integral Block,
         std::copy_constructible Allocator>
constexpr bool
XbitVector<N, Block, Allocator>::invariants() const {
  if (this->begin_ == nullptr) {
    if (this->size_ != 0 || this->cap_ != 0)
      return false;
  } else {
    if (this->cap_ == 0)
      return false;
    if (this->size_ > this->capacity())
      return false;
  }
  return true;
}

}  // namespace biovoltron::detail

namespace biovoltron {

/**
 * @ingroup container
 * Space-efficient container for dibit values, you can 
 * think of it as a specialization of `std::vector` for 
 * type `uint2_t`.
 *
 * The implementation detail is refered to [GCC][GCC] and 
 * [Clang][Clang]'s implementation of [vector<bool>][vector_of_bool] 
 * with adopting new C++20 features such as 
 * [`iterator concepts`][iterator_concept], `operator<=>` etc.
 *
 * The member types and functions definition is just the same as 
 * [`vector<bool>`][vector_of_bool] with following additions:
 * ```cpp
 * block_type; // alias to Block.
 * size_type num_blocks(); // return the number of underlying blocks. 
 * block_type* data();
 * block_type* data() const; // return the begin pointer to the underlying blocks.
 * void flip(); // flip all the dibits/quadbits of the vector.
 * ```
 *
 * Like `vector<bool>`, it can work with all algorithms in 
 * [`<algorithm>`][algorithm] even if [`ranges::sort`][ranges_sort] 
 * which cannot sort `vector<bool>` currently. Note that the `value_type` 
 * of those two containers is `uint8_t` which is not a printable 
 * character, make sure to cast it to `int` before you print:
 * ```cpp
 * std::cout << static_cast<int>(v.front()) << "\n";
 * std::cout << +v.back() << "\n";
 * ```
 *
 * If you want to serialize the containers, you can use 
 * biovoltron::Serializer which proves a simple interface 
 * to serialize binary archive.
 * 
 * Usage
 * ```cpp
 * #include <algorithm>
 * #include <iostream>
 * #include <ranges>
 * #include <biovoltron/container/xbit_vector.hpp>
 * 
 * int main() {
 *   auto v = biovoltron::DibitVector<>{3, 2, 1, 2, 3, 0, 0, 1, 2};
 *   assert(v.size() == 9);
 *   assert(v.num_blocks() == 3);
 * 
 *   auto base_view = std::views::transform([](auto c){ return "ACGT"[c]; });
 *   auto comp_view = std::views::transform([](auto c){ return 0b11u - c; });
 *   std::cout << "Reverse complement of ";
 *   std::ranges::copy(v | base_view, std::ostream_iterator<char>{std::cout, ""});
 *   std::cout << " is ";
 *   for (auto c : v | std::views::reverse | comp_view | base_view)
 *     std::cout << c;
 * }
 * ```
 * 
 * [GCC]: https://github.com/gcc-mirror/gcc/blob/master/libstdc%2B%2B-v3/include/bits/stl_bvector.h
 * [Clang]: https://github.com/llvm-mirror/libcxx/blob/master/include/__bit_reference
 * [vector_of_bool]: https://en.cppreference.com/w/cpp/container/vector_bool
 * [iterator_concept]: https://en.cppreference.com/w/cpp/iterator
 * [algorithm]: https://en.cppreference.com/w/cpp/algorithm
 * [ranges_sort]: https://godbolt.org/z/xb1195
 * [Biovoltron.Serializer]: https://github.com/hewillk/serializer
 * [godbolt]: https://godbolt.org/z/YM5P61
 *
 * @tparam Block Refers to the underlying storage type, 
 * the default type is `uint8_t` which can store 4 dibits
 */
template<std::unsigned_integral Block = std::uint8_t,
         std::copy_constructible Allocator = std::allocator<Block> >
using DibitVector = detail::XbitVector<2, Block, Allocator>;

/**
 * @ingroup container
 * Space-efficient container for quabit values, you can 
 * think of it as a specialization of `std::vector` for 
 * type `uint4_t`.
 *
 * @tparam Block Refers to the underlying storage type, 
 * the default type is `uint8_t` which can store 2 quadbits.
 */
template<std::unsigned_integral Block = std::uint8_t,
         std::copy_constructible Allocator = std::allocator<Block> >
using QuadbitVector = detail::XbitVector<4, Block, Allocator>;

}  // namespace biovoltron
