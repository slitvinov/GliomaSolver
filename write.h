template <typename TWavelets, typename TBlock, typename TLab>
static int write(MRAG::Grid<TWavelets, TBlock> *inputGrid,
                 MRAG::BoundaryInfo *bInfo, const char *path) {
  int np, nt, i, iz, iy, ix, ixx, iyy, izz, j;
  const int shift[8][3] = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0},
      {1, 0, 0}, {1, 0, 1}, {1, 1, 1}, {1, 1, 0},
  };
  int64_t verts[1 + 8];
  vector<MRAG::BlockInfo> vInfo = inputGrid->getBlocksInfo();
  float x[3], y[3];
  FILE *file;
  np = 0;
  file = fopen("a.xyz.raw", "w");
  for (i = 0; i < vInfo.size(); i++) {
    for (iz = 0; iz < TBlock::sizeZ + 1; iz++)
      for (iy = 0; iy < TBlock::sizeY + 1; iy++)
        for (ix = 0; ix < TBlock::sizeX + 1; ix++) {
          vInfo[i].pos(x, ix, iy, iz);
          y[0] = x[0] - TWavelets::CenteringOffset * vInfo[i].h[0];
          y[1] = x[1] - TWavelets::CenteringOffset * vInfo[i].h[1];
          y[2] = x[2] - TWavelets::CenteringOffset * vInfo[i].h[2];
          if (fwrite(y, sizeof y, 1, file) != 1) {
            fprintf(stderr, "%s:%d: error: fail to write\n", __FILE__,
                    __LINE__);
            return 1;
          }
          np++;
        }
  }
  fclose(file);

  file = fopen("a.topo.raw", "w");
  nt = 0;
  verts[0] = 9;
  for (i = 0; i < vInfo.size(); i++) {
    for (iz = 0; iz < TBlock::sizeZ; iz++)
      for (iy = 0; iy < TBlock::sizeY; iy++)
        for (ix = 0; ix < TBlock::sizeX; ix++) {
          for (j = 0; j < 8; j++) {
            ixx = ix + shift[j][0];
            iyy = iy + shift[j][1];
            izz = iz + shift[j][2];
            verts[j + 1] = i * (TBlock::sizeX + 1) * (TBlock::sizeY + 1) *
                               (TBlock::sizeZ + 1) +
                           izz * (TBlock::sizeX + 1) * (TBlock::sizeY + 1) +
                           iyy * (TBlock::sizeX + 1) + ixx;
          }
          if (fwrite(verts, sizeof verts, 1, file) != 1) {
            fprintf(stderr, "%s:%d: error: fail to write\n", __FILE__,
                    __LINE__);
            return 1;
          }
          nt += 1;
        }
  }
  fclose(file);

  /*
  for(int ichannel = 0; ichannel<nChannels; ++ichannel)
  {
          const bool bIsCellCentered = (TWavelets::CenteringOffset>0.0f);
          char channelName[256];
          sprintf(channelName,"channel%d",ichannel);
          vtkFloatArray * fa = vtkFloatArray::New();
          fa->SetNumberOfComponents(1);
          fa->SetName(channelName);
          if( bIsCellCentered )
                  fa->SetNumberOfTuples(totalNumberOfCells);
          else
                  fa->SetNumberOfTuples(totalNumberOfPoints);

          float * vpts = (float*)fa->GetVoidPointer(0);

          if( bIsCellCentered )
          {
                  for(int i=0; i<vInfo.size(); i++)
                  {
                          BlockInfo& info = vInfo[i];
                          TBlock& block =
  inputGrid.getBlockCollection()[info.blockID];

                          int icount = 0;
                          Real h = info.h[0];
                          for(int iz=0; iz<TBlock::sizeZ; iz++)
                                  for(int iy=0; iy<TBlock::sizeY; iy++)
                                          for(int ix=0; ix<TBlock::sizeX; ix++)
                                          {
                                                  const float rValue = (
  block(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ; *vpts++ = rValue;
                                                  ++icount;
                                          }
                  }

                  uGrid->GetCellData()->AddArray(fa);

          }
          else
          {
                  TLab lab;

                  int steStart[3] ={0,0,0};
                  int steEnd[3] = {2,2,2};

                  lab.prepare(inputGrid.getBlockCollection(), bInfo ,steStart,
  steEnd);

                  for(int i=0; i<vInfo.size(); i++)
                  {
                          BlockInfo& info = vInfo[i];
                          lab.load(info);
                          Real h = info.h[0];
                          int icount = 0;
                          for(int iz=0; iz<TBlock::sizeZ+1; iz++)
                                  for(int iy=0; iy<TBlock::sizeY+1; iy++)
                                          for(int ix=0; ix<TBlock::sizeX+1;
  ix++)
                                          {
                                                  const float rValue = (
  lab(ix,iy,iz).giveMe(ichannel + iChannelStart, h) ) ; *vpts++ = rValue;
                                                  ++icount;
                                          }
                  }
                  uGrid->GetPointData()->AddArray(fa);
          }
          fa->Delete();
  }

  vtkXMLUnstructuredGridWriter * writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename.c_str());
  writer->SetInput(uGrid);
  writer->Write();

  uGrid->Delete();
  writer->Delete();
  */

  file = fopen("a.xdmf2", "w");
  fprintf(file,
          "<Xdmf>\n"
          "  <Domain>\n"
          "    <Grid>\n"
          "      <Geometry>\n"
          "        <DataItem\n"
          "            Dimensions=\"%d 3\"\n"
          "            Format=\"Binary\">\n"
          "          %s\n"
          "        </DataItem>\n"
          "      </Geometry>\n"
          "      <Topology\n"
          "	  Dimensions=\"%d\"\n"
          "	  Type=\"Mixed\">\n"
          "	<DataItem\n"
          "	    DataType=\"Int\"\n"
          "	    Dimensions=\"%d\"\n"
          "	    Format=\"Binary\"\n"
          "	    Precision=\"8\">\n"
          "	  %s\n"
          "	</DataItem>\n"
          "      </Topology>\n"
          "    </Grid>\n"
          "  </Domain>\n"
          "</Xdmf>\n",
          np, "a.xyz.raw", nt, nt * (8 + 1), "a.topo.raw");
  fclose(file);
  return 0;
}
