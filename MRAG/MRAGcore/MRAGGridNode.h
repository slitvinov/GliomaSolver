namespace MRAG {
struct GridNode {
  bool isEmpty;
  GridNode *parent;
  int blockID;
  int index[3];
  int level;
  bool shouldBeCompressed;
  bool shouldBeRefined;

  GridNode(bool isEmpty_, GridNode *parent_, int blockID_, int idX, int idY,
           int idZ, int level_)
      : isEmpty(isEmpty_), parent(parent_), blockID(blockID_), level(level_),
        shouldBeCompressed(false), shouldBeRefined(false) {
    index[0] = idX;
    index[1] = idY;
    index[2] = idZ;
  }
};
} // namespace MRAG
