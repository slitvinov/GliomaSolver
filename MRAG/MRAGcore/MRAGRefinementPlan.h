namespace MRAG {
struct RefinementPlanNode {
  int index[3];
  int relative_index[3];
  int level;
};

struct RefinementReport {
  RefinementReport() : blockID(0), childrenIDs() {}
  int blockID;
  std::vector<int> childrenIDs;
};

struct SingleRefinement {
  // int blockID;
  BlockInfo block_info_source;
  std::vector<RefinementPlanNode *> children;
  // std::vector<int> vBlockIDsToRemove;

  SingleRefinement() : block_info_source(), children() {}

  SingleRefinement(SingleRefinement &c) : block_info_source(), children() {
    abort();
  }

  SingleRefinement(const SingleRefinement &c)
      : block_info_source(), children() {
    abort();
  }

  RefinementPlanNode &createEntry() {
    RefinementPlanNode *result = new RefinementPlanNode;

    children.push_back(result);

    return *result;
  }

  ~SingleRefinement() {
    for (int i = 0; i < children.size(); i++)
      delete children[i];
  }
};

struct RefinementPlan {
  std::vector<SingleRefinement *> refinements;

  int nNewBlocks;

  SingleRefinement *createEntry() {
    SingleRefinement *c = new SingleRefinement;

    refinements.push_back(c);

    return c;
  }

  RefinementPlan() : refinements(), nNewBlocks(0) {}

  ~RefinementPlan() {
    for (int i = 0; i < refinements.size(); i++)
      delete refinements[i];
  }
};

struct RefinementResult {
private:
  bool bFailed;

public:
  int nChildrenBlocks;
  int nCollapsedParentBlocks;

  bool hasFailed() const { return bFailed; }

  RefinementResult(const bool bFailed_ = false)
    : bFailed(bFailed_), nChildrenBlocks(0), nCollapsedParentBlocks(0) {}

  RefinementResult operator+=(const RefinementResult &r) {
    nChildrenBlocks += r.nChildrenBlocks;
    nCollapsedParentBlocks += r.nCollapsedParentBlocks;
    bFailed = bFailed || r.bFailed;

    return *this;
  }
};

struct NodeToRefine {
  const GridNode *node;
  int refinementID;
  NodeToRefine(const GridNode *n, int ID) : node(n), refinementID(ID) {}
};
} // namespace MRAG
