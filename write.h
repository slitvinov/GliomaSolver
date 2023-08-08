template <typename TWavelets, typename TBlock, typename TLab>
static int write(MRAG::Grid<TWavelets, TBlock> *inputGrid,
                 MRAG::BoundaryInfo *bInfo, const char *path) {
  int np, nc, i, iz, iy, ix, ixx, iyy, izz, j;
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
  nc = 0;
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
          nc++;
        }
  }
  fclose(file);
  assert(TWavelets::CenteringOffset == 0);
  float value;
  TLab lab;
  int steStart[3] = {0, 0, 0};
  int steEnd[3] = {2, 2, 2};
  lab.prepare(inputGrid->getBlockCollection(), *bInfo, steStart, steEnd);
  file = fopen("a.attr.raw", "w");
  for (i = 0; i < vInfo.size(); i++) {
    MRAG::BlockInfo &info = vInfo[i];
    lab.load(info);
    for (iz = 0; iz < TBlock::sizeZ + 1; iz++)
      for (iy = 0; iy < TBlock::sizeY + 1; iy++)
        for (ix = 0; ix < TBlock::sizeX + 1; ix++) {
          value = lab(ix, iy, iz).phi;
          if (fwrite(&value, sizeof value, 1, file) != 1) {
            fprintf(stderr, "%s:%d: error: fail to write\n", __FILE__,
                    __LINE__);
            return 1;
          }
        }
  }
  fclose(file);

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
          "	  <DataItem\n"
          "	      DataType=\"Int\"\n"
          "	      Dimensions=\"%d\"\n"
          "	      Format=\"Binary\"\n"
          "	      Precision=\"8\">\n"
          "	    %s\n"
          "	  </DataItem>\n"
          "      </Topology>\n",
          np, "a.xyz.raw", nc, nc * (8 + 1), "a.topo.raw");

  fprintf(file,
          "      <Attribute\n"
          "          Name=\"%s\">\n"
          "        <DataItem\n"
          "	       Dimensions=\"%d\"\n"
          "	       Format=\"Binary\">\n"
          "	    %s\n"
          "	   </DataItem>\n"
          "      </Attribute>\n",
          "phi", np, "a.attr.raw");
  fprintf(file, "    </Grid>\n"
                "  </Domain>\n"
                "</Xdmf>\n");
  fclose(file);
  return 0;
}
