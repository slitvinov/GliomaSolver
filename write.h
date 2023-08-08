template <typename TWavelets, typename TBlock, typename TLab>
static int write(MRAG::Grid<TWavelets, TBlock> *inputGrid, MRAG::BoundaryInfo *bInfo, const char *path) {
  int np, i, iz, iy, ix;
  vector<MRAG::BlockInfo> vInfo = inputGrid->getBlocksInfo();
  float x[3], y[3];
  FILE *file;
  np = 0;
  file = fopen("a.raw", "w");
  for (i = 0; i < vInfo.size(); i++) {
    for (iz = 0; iz < TBlock::sizeZ + 1; iz++)
      for (iy = 0; iy < TBlock::sizeY + 1; iy++)
	for (ix = 0; ix < TBlock::sizeX + 1; ix++) {
	  vInfo[i].pos(x, ix, iy, iz);
	  y[0] = x[0] - TWavelets::CenteringOffset * vInfo[i].h[0];
	  y[1] = x[1] - TWavelets::CenteringOffset * vInfo[i].h[1];
	  y[2] = x[2] - TWavelets::CenteringOffset * vInfo[i].h[2];
	  if (fwrite(y, sizeof y, 1, file) != 1) {
	    fprintf(stderr, "%s:%d: error: fail to write\n", __FILE__, __LINE__);
	    return 1;
	  }
	  np++;
	}
  }
  fclose(file);

  /*
  vtkUnstructuredGrid * uGrid = vtkUnstructuredGrid::New();
  uGrid->SetPoints(points);
  points->Delete();
  uGrid->Allocate(totalNumberOfCells);
  vtkIdType verts[verticesPerCell];

  unsigned int counter = 0;
  for(int i=0; i<vInfo.size(); i++)
  {
	  for(int iz=0; iz<TBlock::sizeZ; iz++)
		  for(int iy=0; iy<TBlock::sizeY; iy++)
			  for(int ix=0; ix<TBlock::sizeX; ix++)
			  {
				  for(int j = 0; j < verticesPerCell; j++)
				  {
					  int shift[3] = { j&1, (j>>1)&1,
  (j>>2)&1 }; int ixx = ix + shift[0]; int iyy = iy + shift[1]; int izz = iz +
  shift[2]; int pointIndex =
  counter*(TBlock::sizeX+1)*(TBlock::sizeY+1)*(TBlock::sizeZ+1) +
  izz*(TBlock::sizeX+1)*(TBlock::sizeY+1) + iyy*(TBlock::sizeX+1) + ixx;
					  verts[j] = pointIndex;
				  }
				  uGrid->InsertNextCell(VTK_VOXEL,
  verticesPerCell, verts);
			  }
	  counter += 1;
  }

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
	  "    </Grid>\n"
	  "  </Domain>\n"
	  "</Xdmf>\n",
	  np, "a.raw");
  fclose(file);
  return 0;
}
