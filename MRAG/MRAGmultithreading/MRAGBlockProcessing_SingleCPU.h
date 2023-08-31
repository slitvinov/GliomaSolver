namespace MRAG {
namespace Multithreading {
template <typename BlockType, typename ProcessingMT>
class BlockProcessingMT_Simple {
  ProcessingMT &processing;
  std::vector<BlockInfo> &vInfo;
  BlockType **ptrBlocks;

public:
  BlockProcessingMT_Simple(std::vector<BlockInfo> &vInfo_,
                           ProcessingMT &processing_, BlockType **ptrs)
      : processing(processing_), vInfo(vInfo_), ptrBlocks(ptrs) {}

  template <typename BlockedRange>
  void operator()(const BlockedRange &r) const {
    typedef BlockType B;
    typedef typename B::ElementType E;

    const int nBlocks = r.end() - r.begin();
    BlockInfo *v = new BlockInfo[nBlocks];

    for (int i = 0; i < nBlocks; i++)
      v[i] = vInfo[i + r.begin()];

    // copy constructor required for ProcessingMT
    ProcessingMT p = processing;

    for (int iB = 0; iB < nBlocks; iB++) {
      // operator()(const BlockInfo&, BlockType&) required for ProcessingMT
      p(v[iB], *ptrBlocks[iB + r.begin()]);
    }

    delete[] v;
  }
};

/**
 * Functor to actually perform the operations on the blocks.
 * See MRAG::Multithreading::DummyBlockFunctor for a sample ProcessingMT type.
 */
template <typename BlockType, template <typename BB> class Lab,
          typename Collection, typename ProcessingMT>
class BlockProcessingMT {
  std::vector<BlockInfo> &vInfo;
  Collection &collection;
  BoundaryInfo &boundaryInfo;
  ProcessingMT &processing;
  BlockType **ptrBlocks;

public:
  BlockProcessingMT(std::vector<BlockInfo> &vInfo_, Collection &collection_,
                    BoundaryInfo &boundaryInfo_, ProcessingMT &processing_,
                    BlockType **ptrs)
      : vInfo(vInfo_), collection(collection_), boundaryInfo(boundaryInfo_),
        processing(processing_), ptrBlocks(ptrs) {}

  template <typename BlockedRange>
  void operator()(const BlockedRange &r) const {
    typedef BlockType B;
    typedef typename B::ElementType E;

    const int nBlocks = r.end() - r.begin();
    BlockInfo *v = new BlockInfo[nBlocks];
    for (int i = 0; i < nBlocks; i++)
      v[i] = vInfo[i + r.begin()];

    Lab<B> lab;
    // copy constructor required for ProcessingMT
    ProcessingMT p = processing;
    // stencil_start and stencil_end required for ProcessingMT
    lab.prepare(collection, boundaryInfo, processing.stencil_start,
                processing.stencil_end);

    for (int iB = 0; iB < nBlocks; iB++) {
      BlockInfo &info = v[iB];
      B &b = *ptrBlocks[r.begin() + iB];
      lab.load(info);
      // operator()(LabType&, const BlockInfo&, BlockType&) required for
      // ProcessingMT
      p(lab, info, b);
    }

    delete[] v;
  }
};

/**
 * Process blocks on a single core.
 * @see BlockProcessing_SingleCPU::process
 */
template <typename BlockType> class BlockProcessing_SingleCPU {
  template <typename Collection>
  static BlockType **createBlockPointers(std::vector<BlockInfo> &vInfo,
                                         Collection &collection) {
    BlockType **result = new BlockType *[vInfo.size()];

    for (int i = 0; i < vInfo.size(); i++)
      result[i] = &collection.lock(vInfo[i].blockID);

    return result;
  }

  template <typename Collection>
  static void destroyBlockPointers(BlockType **&ptrs,
                                   std::vector<BlockInfo> &vInfo,
                                   Collection &collection) {
    for (int i = 0; i < vInfo.size(); i++)
      collection.release(vInfo[i].blockID);

    delete[] ptrs;
    ptrs = NULL;
  }

public:
  /**
   * Process blocks one after the other.
   * @param vInfo         Info of all the blocks (to be processed) in the grid
   *                      (e.g. result of Grid::getBlocksInfo()).
   * @param c             Collection of all the blocks (to be processed) in the
   * grid (e.g. result of Grid::getBlockCollection()).
   * @param p             Functor processing the block.
   *                      See MRAG::Multithreading::DummySimpleBlockFunctor for
   * details.
   * @param dummy         Optional and never used (just to match signature of
   * tbb-versions).
   */
  template <typename Processing, typename Collection>
  static void process(std::vector<BlockInfo> &vInfo, Collection &c,
                      Processing &p, int dummy = -1) {
    BlockType **ptrs = createBlockPointers(vInfo, c);

    BlockProcessingMT_Simple<BlockType, Processing> bps(vInfo, p, ptrs);
    bps(SimpleInterval(0, vInfo.size()));

    destroyBlockPointers(ptrs, vInfo, c);
  }

  /**
   * Process blocks one after the other.
   * @param vInfo         Info of all the blocks (to be processed) in the grid
   *                      (e.g. result of Grid::getBlocksInfo()).
   * @param c             Collection of all the blocks (to be processed) in the
   * grid (e.g. result of Grid::getBlockCollection()).
   * @param b             Info on the boundaries of the grid (e.g. result of
   * Grid::getBoundaryInfo())
   * @param p             Functor processing the block.
   *                      See MRAG::Multithreading::DummyBlockFunctor for
   * details.
   * @param dummy         Optional and never used (just to match signature of
   * tbb-versions).
   */
  template <template <typename Btype> class Lab, typename Processing,
            typename Collection>
  static void process(std::vector<BlockInfo> &vInfo, Collection &c,
                      BoundaryInfo &b, Processing &p, int dummy = -1) {
    BlockType **ptrs = createBlockPointers(vInfo, c);

    BlockProcessingMT<BlockType, Lab, Collection, Processing> bpg(vInfo, c, b,
                                                                  p, ptrs);
    bpg(SimpleInterval(0, vInfo.size()));

    destroyBlockPointers(ptrs, vInfo, c);
  }

  /**
   * Alias for process<BlockLab>(vInfo, c, b, p).
   * @see BlockProcessing_SingleCPU::process()
   */
  template <typename Processing, typename Collection>
  void pipeline_process(std::vector<BlockInfo> &vInfo, Collection &c,
                        BoundaryInfo &b, Processing &p) {
    process<BlockLab>(vInfo, c, b, p);
  }

  template <typename B> friend class BlockProcessing_TBB;
};
} /* namespace Multithreading */
} /* namespace MRAG */
