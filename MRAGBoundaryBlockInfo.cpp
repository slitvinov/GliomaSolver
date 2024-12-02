#include <algorithm>
#include <assert.h>
#include <map>
#include <math.h>
#include <string.h>
#include <queue>
#include <stack>
#include <stdio.h>
#include "MRAG/MRAGcore/MRAGCommon.h"
#include "MRAG/MRAGcore/MRAGBitStream.h"
#include "MRAG/MRAGcore/MRAGEncoder.h"
#include "MRAG/MRAGcore/MRAGHuffmanEncoder.h"
#include "MRAG/MRAGcore/MRAGBoundaryBlockInfo.h"
namespace MRAG {
void BoundaryInfoBlock::_compress() {
  // 1. encode number of instruction per ghosts
  // 2. encode points
  // 3. encode weights indices

  int nItems = 0;
  int nMaxInstrPerGhost = 0;

  // 1.
  {
    // A. iterate over the ghosts and retrieve the number of instruction items
    // for each ghost (and fill nItems) B. compress that

    std::vector<unsigned short> vStream(ghosts.size());

    std::vector<ReconstructionInfo>::iterator itSource = ghosts.begin();
    const std::vector<unsigned short>::iterator itEnd = vStream.end();

    // A.
    for (std::vector<unsigned short>::iterator it = vStream.begin(); it != itEnd;
         it++, itSource++) {
      const int s = itSource->size();
      *it = s;
      nItems += s;
      nMaxInstrPerGhost = std::max(nMaxInstrPerGhost, s);
    }

    // B.
    encodedInstructionSizes.encode(vStream, nMaxInstrPerGhost + 1);
  }

  std::vector<int> vIP2newIP(indexPool.size());

  // 2.
  {
    // A. group points by blockID
    // B. count/store how many points there are for each blockID in the
    // indexpool C. iterating on each blockID, stream the 3Dindexed point_index
    // D. encode that

    std::map<int, std::vector<unsigned short>> mapBlockID2index;

    // A.
    {
      const std::vector<PointIndex>::const_iterator itEnd = indexPool.end();
      int c = 0;
      for (std::vector<PointIndex>::const_iterator it = indexPool.begin();
           it != itEnd; it++, c++) {
        std::vector<unsigned short> &v = mapBlockID2index[it->blockID];
        vIP2newIP[c] = v.size();
        v.push_back(it->index);
      }

      std::map<int, int> mapCumulativeOffset;

      {
        int s = 0;
        const std::map<int, std::vector<unsigned short>>::iterator itEnd =
            mapBlockID2index.end();
        for (std::map<int, std::vector<unsigned short>>::iterator it =
                 mapBlockID2index.begin();
             it != itEnd; it++) {
          mapCumulativeOffset[it->first] = s;
          s += it->second.size();
        }
      }

      {
        int c = 0;
        for (std::vector<PointIndex>::const_iterator it = indexPool.begin();
             it != itEnd; it++, c++)
          vIP2newIP[c] += mapCumulativeOffset[it->blockID];
      }
    }

    // B.
    {
      assert(vBlockID_Points.size() == 0);
      vBlockID_Points.resize(mapBlockID2index.size());
      std::map<int, std::vector<unsigned short>>::const_iterator itSource =
          mapBlockID2index.begin();
      const std::vector<std::pair<int, int>>::iterator itEnd = vBlockID_Points.end();
      for (std::vector<std::pair<int, int>>::iterator it = vBlockID_Points.begin();
           it != itEnd; it++, itSource++) {
        it->first = itSource->first;
        it->second = itSource->second.size();
      }
    }

    std::vector<unsigned char> vStream(indexPool.size() * 3);
    int maxIndexValue = 0;

    // C.
    {
      std::vector<unsigned char>::iterator itDest = vStream.begin();

      const std::map<int, std::vector<unsigned short>>::const_iterator itVEnd =
          mapBlockID2index.end();
      for (std::map<int, std::vector<unsigned short>>::const_iterator itV =
               mapBlockID2index.begin();
           itV != itVEnd; itV++) {
        std::vector<unsigned short>::const_iterator itEnd = itV->second.end();
        for (std::vector<unsigned short>::const_iterator itS = itV->second.begin();
             itS != itEnd; itS++) {
          const int index = *itS;
          const unsigned char i3D[3] = {
              (unsigned char)((index % block_size[0])),
              (unsigned char)((index / block_size[0]) % block_size[1]),
              (unsigned char)(index / (block_size[0] * block_size[1]))};

          assert(i3D[2] < block_size[2]);

          for (int i = 0; i < 3; i++)
            maxIndexValue = std::max(maxIndexValue, (int)i3D[i]);

          for (int i = 0; i < 3; i++, itDest++)
            *itDest = i3D[i];
        }
      }
    }

    // D.
    vBlockID_encodedPointIndices3D.encode(vStream, maxIndexValue + 1);
  }

  // 3.
  {
    // A. create a stream of separated, meaningless indices from
    // weights_index[3] B. encode that C. create a stream of indices (referring
    // to the index pool) D. encode that

    // A.-B.
    {
      std::vector<unsigned char> streamW(nItems * 3);
      std::vector<unsigned char>::iterator itDest = streamW.begin();

      // A.
      const std::vector<ReconstructionInfo>::const_iterator itGEnd = ghosts.end();
      for (std::vector<ReconstructionInfo>::const_iterator itG = ghosts.begin();
           itG != itGEnd; itG++) {
        const std::vector<IndexWP>::const_iterator itItemEnd = itG->end();

        for (std::vector<IndexWP>::const_iterator itItem = itG->begin();
             itItem != itItemEnd; itItem++)
          for (int i = 0; i < 3; i++, itDest++)
            *itDest = itItem->weights_index[i];
      }

      // B.
      encodedInstructionItemsWs.encode(streamW, weightsPool.size() + 1);
    }

    // C.-D.
    {
      std::vector<unsigned short> streamI(nItems);
      std::vector<unsigned short>::iterator itDest = streamI.begin();

      // C.
      const std::vector<ReconstructionInfo>::const_iterator itGEnd = ghosts.end();
      for (std::vector<ReconstructionInfo>::const_iterator itG = ghosts.begin();
           itG != itGEnd; itG++) {
        const std::vector<IndexWP>::const_iterator itItemEnd = itG->end();

        for (std::vector<IndexWP>::const_iterator itItem = itG->begin();
             itItem != itItemEnd; itItem++, itDest++)
          *itDest = vIP2newIP[itItem->point_index];
      }

      // D.
      encodedInstructionItemsPts.encode(streamI, indexPool.size());
    }
  }
}

void BoundaryInfoBlock::_decompress() {
  // 1. build empty ghost vectors from encodedInstructionSizes
  // 2. build ghost instruction items from encodedInstructionItemsWs and
  // encodedInstructionItemsPts
  // 3. build indexPool from vBlockID_encodedPointIndices3D and vBlockID_Points

  // 1.
  {
    assert(ghosts.size() == 0);

    std::vector<unsigned short> vStream;
    encodedInstructionSizes.decode(vStream);

    ghosts.resize(vStream.size());

    std::vector<ReconstructionInfo>::iterator itDest = ghosts.begin();
    const std::vector<unsigned short>::const_iterator itStreamEnd = vStream.end();
    for (std::vector<unsigned short>::const_iterator itStream = vStream.begin();
         itStream != itStreamEnd; itStream++, itDest++)
      itDest->resize(*itStream);
  }

  // 2.
  {// A. decode the indices for the weights pool
   // B. decode the indices for the index pool

   // A.
   {std::vector<unsigned char> vIW;
  encodedInstructionItemsWs.decode(vIW);

  std::vector<unsigned char>::iterator itSource = vIW.begin();
  const std::vector<ReconstructionInfo>::const_iterator itGEnd = ghosts.end();
  for (std::vector<ReconstructionInfo>::iterator itG = ghosts.begin(); itG != itGEnd;
       itG++) {
    const std::vector<IndexWP>::iterator itItemEnd = itG->end();

    for (std::vector<IndexWP>::iterator itItem = itG->begin(); itItem != itItemEnd;
         itItem++)
      for (int i = 0; i < 3; i++, itSource++)
        itItem->weights_index[i] = *itSource;
  }
}

// B.
{
  std::vector<unsigned short> vIP;
  encodedInstructionItemsPts.decode(vIP);

  std::vector<unsigned short>::iterator itSource = vIP.begin();
  const std::vector<ReconstructionInfo>::const_iterator itGEnd = ghosts.end();
  for (std::vector<ReconstructionInfo>::iterator itG = ghosts.begin(); itG != itGEnd;
       itG++) {
    const std::vector<IndexWP>::iterator itItemEnd = itG->end();

    for (std::vector<IndexWP>::iterator itItem = itG->begin(); itItem != itItemEnd;
         itItem++, itSource++)
      itItem->point_index = *itSource;
  }
}
} // namespace MRAG

// 3.
{
  assert(indexPool.size() == 0);

  std::vector<unsigned char> vStream;

  vBlockID_encodedPointIndices3D.decode(vStream);

  assert(vStream.size() % 3 == 0);
  indexPool.resize(vStream.size() / 3);

  std::vector<PointIndex>::iterator itDest = indexPool.begin();
  std::vector<unsigned char>::iterator itSource = vStream.begin();

  const int row_size = block_size[0];
  const int slice_size = block_size[0] * block_size[1];

  const std::vector<std::pair<int, int>>::const_iterator itVEnd = vBlockID_Points.end();
  for (std::vector<std::pair<int, int>>::const_iterator itV = vBlockID_Points.begin();
       itV != itVEnd; itV++) {
    const int blockID = itV->first;
    const int nPoints = itV->second;

    for (int i = 0; i < nPoints; i++, itDest++) {
      unsigned char i3D[3] = {0, 0, 0};

      for (int i = 0; i < 3; i++, itSource++)
        i3D[i] = *itSource;

      const int index = i3D[0] + i3D[1] * row_size + i3D[2] * slice_size;

      itDest->blockID = blockID;
      itDest->index = index;
    }
  }
}
}

void BoundaryInfoBlock::_discard_decompression() {
  // 1. discard indexPool content
  // 2. discard ghost instruction content

  // 1.
  indexPool.clear();

  // 2.
  const std::vector<ReconstructionInfo>::iterator itGEnd = ghosts.end();
  for (std::vector<ReconstructionInfo>::iterator itG = ghosts.begin(); itG != itGEnd;
       itG++)
    itG->clear();

  ghosts.clear();
}

void *BoundaryInfoBlock::createBBPack() { return NULL; }

float BoundaryInfoBlock::getMemorySize() const {
  int memsize = sizeof(BoundaryInfoBlock);

  memsize += indexPool.size() * sizeof(PointIndex);

  memsize += weightsPool.size() * sizeof(double);

  memsize += dependentBlockIDs.size() * sizeof(int);

  int ghost_size = 0;
  for (int i = 0; i < ghosts.size(); i++)
    ghost_size += ghosts[i].size() * sizeof(IndexWP);

  memsize += ghost_size;

  double compressedDataMB =
      encodedInstructionSizes.getMemorySize() +
      encodedInstructionItemsWs.getMemorySize() +
      encodedInstructionItemsPts.getMemorySize() +
      vBlockID_encodedPointIndices3D.getMemorySize() +
      vBlockID_Points.size() * sizeof(std::pair<int, int>) / 1024. / 1024.;

  // printf("M*** MEMSIIZE UNC = %d MEMSIZE COMP=%d\n", ghost_size,
  // (int)(compressedDataMB*1024.*1024.));
  memsize += compressedDataMB * 1024. * 1024.;

  const bool bVerbose = false;

  if (bVerbose)
    printf("BoundaryInfoBlock MB=%f, indexPool=%f%%, weightsPool=%f%%, "
           "instruction=%f%%\n",
           memsize / (double)(1 << 20),
           indexPool.size() * sizeof(PointIndex) / (double)memsize * 100.0,
           weightsPool.size() * sizeof(double) / (double)memsize * 100.0,
           ghost_size / (double)memsize * 100);

  return memsize / (double)(1 << 20);
}

float BoundaryInfo::getMemorySize() const {
  int memsize = sizeof(BoundaryInfo) + sizeof(int) * boundaryInfoOfBlock.size();

  double sum = 0;

  for (std::map<int, BoundaryInfoBlock *>::const_iterator it =
           boundaryInfoOfBlock.begin();
       it != boundaryInfoOfBlock.end(); it++)
    sum += it->second->getMemorySize();

  return (sum + memsize / (double)(1 << 20));
}

void BoundaryInfo::clear() {
  for (std::map<int, BoundaryInfoBlock *>::iterator it = boundaryInfoOfBlock.begin();
       it != boundaryInfoOfBlock.end(); it++) {
    delete it->second;

    it->second = NULL;
  }
  boundaryInfoOfBlock.clear();
}

void BoundaryInfo::erase(int blockID, bool bCritical) {
  std::map<int, BoundaryInfoBlock *>::iterator it =
      boundaryInfoOfBlock.find(blockID);

  assert(!bCritical || it != boundaryInfoOfBlock.end());

  if (!bCritical && it == boundaryInfoOfBlock.end())
    return;

  delete it->second;

  it->second = NULL;

  boundaryInfoOfBlock.erase(blockID);
}
}
