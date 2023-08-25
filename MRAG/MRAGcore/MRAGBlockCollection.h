namespace MRAG {

/**
 * Collection of blocks (used by MRAG::Grid to collect blocks).
 */
template <typename BlockType_> class BlockCollection {
public:
  static const int nChunkSize =
      1; // 2*(BlockType_::shouldProcessDirectionY?2:1)*(BlockType_::shouldProcessDirectionZ?2:1)
         // - 1;
  typedef BlockType_ BlockType;

  /** Default constructor creating empty collection. */
  BlockCollection();
  /** Destructor making sure, we clean up all blocks. */
  ~BlockCollection();

  /**
   * Create or resize collection of blocks to contain certain number of blocks.
   * @param nNumberOfBlocks   Number of blocks to have in this collection.
   * @return Vector of IDs to access the blocks.
   */
  std::vector<int> create(const int nNumberOfBlocks = 1);
  /** Remove a certain block. */
  void erase(const int ID);
  /** Resets collection to an empty one. */
  void clear();
  /** Access a certain block. */
  virtual BlockType &operator[](const int blockID) const;

  // REMOTE ACCESS INTERFACE
  virtual BlockType &lock(const int blockID) const { return (*this)[blockID]; }

  virtual void release(const int blockID, const bool bWriteBack = true) const {}

  virtual void lock(const std::vector<int> &blockIDs,
                    std::vector<BlockType *> *output = NULL) const {
    if (output != NULL) {
      output->resize(blockIDs.size());

      for (int i = 0; i < blockIDs.size(); i++)
        (*output)[i] = &lock(blockIDs[i]);
    } else
      for (std::vector<int>::const_iterator it = blockIDs.begin();
           it != blockIDs.end(); it++)
        lock(*it);
  }

  void release(const std::vector<int> &blockIDs, bool bWriteBack = true) const {
    for (std::vector<int>::const_iterator it = blockIDs.begin();
         it != blockIDs.end(); it++)
      release(*it);
  }

  virtual float getMemorySize(bool bCountTrashAlso = false) const;

  float getBlockSize() const { return sizeof(BlockType); }

protected:
  struct Chunk {
    BlockType *p;
    int startID;
    int nActives;

    Chunk(BlockType *p_, const int n_, const int startID_)
        : p(p_), startID(startID_), nActives(n_) {}
    ~Chunk() {}

    Chunk(const Chunk &a) : p(a.p), startID(a.startID), nActives(a.nActives) {}

    Chunk &operator=(const Chunk &a) {
      p = (a.p);
      startID = (a.startID);
      nActives = (a.nActives);

      return *this;
    }
  };

  virtual BlockType *_allocate(int nElements) const {
    BlockType *ptr = allocator.allocate(nElements);

    for (int i = 0; i < nElements; i++)
      allocator.construct(ptr + i, BlockType_());

    return ptr;
  }

  virtual void _deallocate(BlockType *ptr, int nElements) const {
    for (int i = 0; i < nElements; i++)
      allocator.destroy(ptr + i);

    allocator.deallocate(ptr, nElements);
  }

  static int _createIDs(int n = 1);

  std::vector<int>
  _allocateBlockInChunk(int &inoutRequestedBlocks, Chunk *chunk,
                        int &inoutAvailableBlocksInCurrentChunk) const;
  std::vector<int> _allocateChunks(int &inoutRequestedBlocks,
                              set<Chunk *> &inoutChunks,
                              std::map<int, BlockType *> &inoutIDToBlockPointers,
                              std::map<int, Chunk *> &inoutBlockIDToChunck,
                              std::vector<Chunk *> &inoutRecycledChunks,
                              Chunk *&outCurrentChunk,
                              int &outAvailableBlocksInCurrentChunk) const;
  void _emptyTrash(int nChunks = -1);
  void _trash();

  std::map<int, BlockType *> m_blockIDToBlockPointers;
  std::map<int, Chunk *> m_blockIDToChunck;
  set<Chunk *> m_setChunks;
  std::vector<Chunk *> m_trash;

  Chunk *m_currentChunk;
  int m_nAvailableBlocksInCurrentChunk;

  mutable std::allocator<BlockType_> allocator;

private:
  // forbidden
  BlockCollection(const BlockCollection &)
      : m_blockIDToBlockPointers(), m_blockIDToChunck(), m_setChunks(),
        m_trash(), m_currentChunk(NULL), m_nAvailableBlocksInCurrentChunk(0) {
    abort();
  }

  BlockCollection &operator=(BlockCollection &) {
    abort();
    return *this;
  }
};

} // namespace MRAG
