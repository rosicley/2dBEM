#include "Problem.h"

Problem::Problem()
{
    char buff[FILENAME_MAX];
    getCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    current_working_dir_ = current_working_dir;
}

Problem::Problem(const int &nIP, const double &collocParam, const double &offsetSourceP)
{
    char buff[FILENAME_MAX];
    getCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    current_working_dir_ = current_working_dir;
    quadrature_ = new Gauss1D(nIP);
    collocParam_ = collocParam;
    if (fabs(offsetSourceP) <= 1.0e-06)
    {
        sourceOut_ = false;
    }
    else
    {
        sourceOut_ = true;
    }
    offsetSourceOut_ = offsetSourceP;
}

Problem::~Problem()
{
}

void Problem::generateMesh(Geometry *geometry, const int &order, const std::string &algorithm, const bool &plotMesh, const bool &deleteFiles, const bool &showInfo)
{
    geometry_ = geometry;
    order_ = order;

    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    std::string file = geometry_->getName();
    std::string geoFile = file + ".geo";
    std::string mshFile = file + ".msh";

    if (rank == 0)
    {
        std::string gmshCode = geometry_->createGmshCode();

        gmshCode += "Mesh.ElementOrder = " + std::to_string(order) + ";\n//\n";
        if (algorithm == "AUTO")
        {
            gmshCode += "Mesh.Algorithm = 2;\n//\n";
        }
        else if (algorithm == "DELAUNAY")
        {
            gmshCode += "Mesh.Algorithm = 5;\n//\n";
        }
        else if (algorithm == "FRONT")
        {
            gmshCode += "Mesh.Algorithm = 6;\n//\n";
        }
        else if (algorithm == "ADAPT")
        {
            gmshCode += "Mesh.Algorithm = 1;\n//\n";
        }
        else if (algorithm == "BAMG")
        {
            gmshCode += "Mesh.Algorithm = 7;\n//\n";
        }
        else
        {
            std::cout << algorithm << " is not supported. Please select another type of algorithm. \n ";
            exit(EXIT_FAILURE);
        }

        gmshCode += "Mesh 1;\n//\n";

        gmshCode += "Mesh.MshFileVersion = 2.2;\n";
        gmshCode += "Save \"" + mshFile + "\";\n";

        std::ofstream file(geoFile);
        file << gmshCode;
        file.close();

        std::string cmd = current_working_dir_ + "/src/mesh/gmsh ";
        if (plotMesh)
        {
            cmd += geoFile; // + " &";
        }
        else
        {
            cmd += geoFile + " -";
        }

        if (!showInfo and !plotMesh)
        {
            cmd += " -v 0";
        }

        system(cmd.c_str());
    }
    MPI_Barrier(PETSC_COMM_WORLD);

    gmshReading(mshFile);

    MPI_Barrier(PETSC_COMM_WORLD);

    if (rank == 0 and deleteFiles)
    {
        std::string aux = geoFile + " " + mshFile;
        system((remove + aux).c_str());
    }
}

void Problem::gmshReading(const std::string &gmshFile)
{
    int nentities, nnode, nelem, trash, auxconec, entitie;
    double coord1, coord2, coord3;

    std::unordered_map<int, std::string> physicalEntities;
    std::string line;
    std::vector<BoundaryElement *> boundAux;

    //opening .msh file
    std::ifstream file(gmshFile);

    int nConec = order_ + 1;
    std::vector<Node *> boundaryConec(nConec);

    jumpLine(4, file);
    file >> nentities;
    std::getline(file, line);
    physicalEntities.reserve(nentities);
    for (int i = 0; i < nentities; i++)
    {
        std::getline(file, line);
        std::vector<std::string> tokens = split(line, " ");
        int index;
        std::istringstream(tokens[1]) >> index;
        physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
    }
    for (std::unordered_map<int, std::string>::const_iterator entities = physicalEntities.begin(); entities != physicalEntities.end(); entities++)
    {
        std::string name = entities->second;
        if (name[0] == 'l')
        {
            lineElements_.insert(std::unordered_map<std::string, std::vector<BoundaryElement *>>::value_type(name, boundAux));
        }
    }

    jumpLine(2, file);
    //reading nodes
    file >> nnode;
    std::getline(file, line);
    for (int inode = 0; inode < nnode; inode++)
    {
        file >> trash >> coord1 >> coord2 >> coord3;
        addNode(inode, {coord1, coord2, coord3});
        std::getline(file, line);
    }

    jumpLine(2, file);
    //reading elements
    file >> nelem;
    std::getline(file, line);
    int cont = 0;
    for (int ielem = 0; ielem < nelem; ielem++)
    {
        file >> trash >> trash >> trash >> entitie >> trash;
        std::string name = physicalEntities[entitie];
        if (name[0] == 'l')
        {
            for (int in = 0; in < nConec; in++)
            {
                file >> auxconec;
                boundaryConec[in] = nodes_[auxconec - 1];
            }
            addBoundaryElement(cont, name, boundaryConec);
            cont++;
        }
        else if (name[0] == 'p' and geometry_->getPoints().count(name) >= 1)
        {
            file >> auxconec;
            geometry_->getPoint(name)->addNodeToPoint(nodes_[auxconec - 1]);
            if (geometry_->verifyDuplicatedPoint(name))
            {
                discontinuousNodes_.insert(nodes_[auxconec - 1]);
            }
        }
        std::getline(file, line);
    }
    file.close();
}

void Problem::addNode(const int &index, const vecDouble &coordinate)
{
    Node *no = new Node(index, coordinate);
    nodes_.push_back(no);

    CollocationPoint *colloc = new CollocationPoint(no);
    collocPoints_.push_back(colloc);
}

void Problem::addBoundaryElement(const int &index, const std::string &name, const std::vector<Node *> &connection)
{
    int node0 = discontinuousNodes_.count(connection[0]);
    int node1 = discontinuousNodes_.count(connection[1]);
    vecDouble coord(3);
    double aux = 2.0 / static_cast<double>(order_) * collocParam_;

    std::vector<CollocationPoint *> collocConec;
    for (Node *n : connection)
    {
        collocConec.push_back(collocPoints_[n->getIndex()]);
    }

    if (node0 == 0 and node1 == 0) //continuous element
    {
        BoundaryElement *bFEM = new BoundaryElement(index, connection, collocConec, quadrature_);
        lineElements_[name].push_back(bFEM);
        elements_.push_back(bFEM);
    }
    else if (node0 == 1 and node1 == 1) //discontinuous element
    {
        BoundaryElement *bFEM = new DiscontBoundaryElement(index, connection, collocConec, quadrature_, "both", collocParam_);
        lineElements_[name].push_back(bFEM);
        elements_.push_back(bFEM);

        bFEM->interpolateGeometricalCoordinate(-1.0 + aux, coord);
        collocPoints_[connection[0]->getIndex()]->setCoordinate(coord);
        bFEM->interpolateGeometricalCoordinate(1.0 - aux, coord);
        collocPoints_[connection[1]->getIndex()]->setCoordinate(coord);
    }
    else if (node0 == 1) //semi-continuous element
    {
        BoundaryElement *bFEM = new DiscontBoundaryElement(index, connection, collocConec, quadrature_, "left", collocParam_);
        lineElements_[name].push_back(bFEM);
        elements_.push_back(bFEM);

        bFEM->interpolateGeometricalCoordinate(-1.0 + aux, coord);
        collocPoints_[connection[0]->getIndex()]->setCoordinate(coord);
    }
    else if (node1 == 1) //semi-continuous element
    {
        BoundaryElement *bFEM = new DiscontBoundaryElement(index, connection, collocConec, quadrature_, "right", collocParam_);
        lineElements_[name].push_back(bFEM);
        elements_.push_back(bFEM);

        bFEM->interpolateGeometricalCoordinate(1.0 - aux, coord);
        collocPoints_[connection[1]->getIndex()]->setCoordinate(coord);
    }
    else
    {
        std::cout << "ERROR...\n";
    }
}

void Problem::createSourcePoints()
{
    for (CollocationPoint *colloc : collocPoints_)
    {
        SourcePoint *source = new SourcePoint(colloc);
        sourcePoints_.push_back(source);
    }

    if (!sourceOut_)
    {
        for (BoundaryElement *el : elements_)
        {
            for (Node *n : el->getCollocationConnection())
            {
                sourcePoints_[n->getIndex()]->addInsideElement(el->getIndex());
            }
        }
    }
    else
    {
        for (BoundaryElement *el : elements_)
        {
            std::vector<Node *> connection = el->getConnection();
            int node0 = discontinuousNodes_.count(connection[0]);
            int node1 = discontinuousNodes_.count(connection[1]);
            vecDouble coord(3);

            vecDouble xsis(order_ + 1);
            xsis[0] = -1.0;
            xsis[1] = 1.0;
            double aux = 2.0 / static_cast<double>(order_);
            for (int i = 0; i < order_ - 1; i++)
            {
                xsis[2 + i] = xsis[0] + static_cast<double>(i + 1) * aux;
            }

            aux = aux * collocParam_;
            if (node0 == 1 and node1 == 1) //discontinuous element
            {
                xsis[0] += aux;
                xsis[1] -= aux;
            }
            else if (node0 == 1)
            {
                xsis[0] += aux;
            }
            else if (node1 == 1)
            {
                xsis[1] -= aux;
            }
            for (int i = 0; i <= order_; i++)
            {
                el->coordOfSourcePoints(xsis[i], offsetSourceOut_, coord);
                sourcePoints_[connection[i]->getIndex()]->setCoordinate(coord);
            }
        }
    }

    for (auto &region : geometry_->getRegions())
    {
        std::vector<BoundaryElement *> elements;
        for (auto &line : region.second)
        {
            for (auto &el : lineElements_[line])
            {
                elements.push_back(el);
            }
        }
        subElements_[region.first] = elements;
    }

    for (auto &region : subElements_)
    {
        std::unordered_set<CollocationPoint *> collocAuxiliar;
        std::vector<SourcePoint *> sourceR;
        for (auto &el : region.second)
        {
            for (auto &colloc : el->getCollocationConnection())
            {
                if (collocAuxiliar.count(colloc) == 0)
                {
                    sourceR.push_back(sourcePoints_[colloc->getIndex()]);
                    collocAuxiliar.insert(colloc);
                }
            }
        }
        subSourcePoints_[region.first] = sourceR;
    }
}

void Problem::exportToParaviewGeometricMesh(const int &index)
{
    std::stringstream nameFile;
    nameFile << "geometricMesh" << index << ".vtu";
    std::ofstream file(nameFile.str());

    int nElement = 0;
    for (auto &el : lineElements_)
    {
        nElement += el.second.size();
    }

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << nElement
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";

    for (Node *n : nodes_)
    {
        file << n->getCoordinate(0) << " " << n->getCoordinate(1) << " " << n->getCoordinate(2) << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";
    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            for (Node *no : el->getConnection())
            {
                file << no->getIndex() << " ";
            }
            file << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            aux += el->getConnection().size();
            file << aux << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 68 << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";
    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Process\" format=\"ascii\">" << std::endl;

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 0 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"LineNumber\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        std::string ln = line.first.substr(1, line.first.size() - 1);
        for (BoundaryElement *el : line.second)
        {
            file << ln << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"ElementIndex\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            file << el->getIndex() << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

void Problem::exportToParaviewCollocationMesh_Potential(const int &index)
{
    std::stringstream nameFile;
    nameFile << "collocationMesh" << index << ".vtu";
    std::ofstream file(nameFile.str());

    int nElement = 0;
    for (auto &el : lineElements_)
    {
        nElement += el.second.size();
    }

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << collocPoints_.size() + internalPoints_.size()
         << "\"  NumberOfCells=\"" << nElement + internalPoints_.size()
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";

    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getCoordinate(0) << " " << n->getCoordinate(1) << " " << n->getCoordinate(2) << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << n->getCoordinate(0) << " " << n->getCoordinate(1) << " " << n->getCoordinate(2) << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";
    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            for (Node *no : el->getConnection())
            {
                file << no->getIndex() << " ";
            }
            file << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << collocPoints_.size() + n->getIndex() << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            aux += el->getConnection().size();
            file << aux << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << ++aux << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 68 << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << 1 << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Potential\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getPotential() << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << internalPotential_[n->getIndex()] << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Flux\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getFlux() << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << sqrt(internalFlux_[n->getIndex()][0] * internalFlux_[n->getIndex()][0] + internalFlux_[n->getIndex()][1] * internalFlux_[n->getIndex()][1]) << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"InternalFlux\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << 0 << " " << 0 << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << internalFlux_[n->getIndex()][0] << " " << internalFlux_[n->getIndex()][1] << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </PointData>"
         << "\n";

    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Process\" format=\"ascii\">" << std::endl;

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 0 << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << 0 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"LineNumber\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        std::string ln = line.first.substr(1, line.first.size() - 1);
        for (BoundaryElement *el : line.second)
        {
            file << ln << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << -1 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"ElementIndex\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            file << el->getIndex() << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << -1 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

void Problem::exportToParaviewSourcePoints(const int &index)
{
    std::stringstream nameFile;
    nameFile << "sourcePoints" << index << ".vtu";
    std::ofstream file(nameFile.str());

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << sourcePoints_.size()
         << "\"  NumberOfCells=\"" << sourcePoints_.size()
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";

    for (auto &sub : subSourcePoints_)
    {
        for (auto &n : sub.second)
        {
            file << n->getCoordinate(0) << " " << n->getCoordinate(1) << " " << n->getCoordinate(2) << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";
    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";
    int cont = 0;
    for (auto &sub : subSourcePoints_)
    {
        for (auto &n : sub.second)
        {
            file << cont++ << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";

    cont = 0;
    for (auto &sub : subSourcePoints_)
    {
        for (auto &n : sub.second)
        {
            file << ++cont << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (auto &sub : subSourcePoints_)
    {
        for (auto &n : sub.second)
        {
            file << 1 << "\n";
        }
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"SubRegion\" format=\"ascii\">"
         << "\n";
    for (auto &sub : subSourcePoints_)
    {
        int index = std::stoi(sub.first.substr(1, sub.first.size() - 0));
        for (auto &n : sub.second)
        {
            file << index << "\n";
        }
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </PointData>"
         << "\n";
    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

int Problem::solvePotentialProblem()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    applyPotentialBoundaryConditions();
    createSourcePoints();

    MPI_Barrier(PETSC_COMM_WORLD);

    int ndegree = collocPoints_.size();

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          ndegree, ndegree, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, ndegree);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &All);
    CHKERRQ(ierr);

    matDouble matGl, matHl;
    if (rank == 0)
    {
        for (BoundaryElement *el : elements_)
        {
            el->potentialContribution(sourcePoints_, matGl, matHl);
            std::vector<CollocationPoint *> conec = el->getCollocationConnection();

            for (int in = 0, nNodes = conec.size(); in < nNodes; in++)
            {
                int indexCP = conec[in]->getIndex();
                if (collocCondition_[indexCP] == 0) ///CONHEÇO O FLUXO NO PONTO
                {
                    for (int isp = 0; isp < ndegree; isp++)
                    {
                        ierr = MatSetValues(A, 1, &isp, 1, &indexCP, &matHl[isp][in], ADD_VALUES); //A é a matriz H e C é a matriz G
                        ierr = MatSetValues(C, 1, &isp, 1, &indexCP, &matGl[isp][in], ADD_VALUES);
                    }
                }
                else ///CONHEÇO O POTENCIAL NO PONTO
                {
                    for (int isp = 0; isp < ndegree; isp++)
                    {
                        double Haux = -matHl[isp][in];
                        double Gaux = -matGl[isp][in];

                        ierr = MatSetValues(A, 1, &isp, 1, &indexCP, &Gaux, ADD_VALUES);
                        ierr = MatSetValues(C, 1, &isp, 1, &indexCP, &Haux, ADD_VALUES);
                    }
                }
            }
        }
    }

    if (!sourceOut_ and rank == 0)
    {
        std::cout << "PONTO FONTE DENTRO\n";
        double ccc = 0.5;
        for (int isp = 0; isp < ndegree; isp++)
        {
            if (collocCondition_[isp] == 0) ///CONHEÇO O FLUXO NO PONTO
            {
                ierr = MatSetValues(A, 1, &isp, 1, &isp, &ccc, ADD_VALUES);
            }
            else ///CONHEÇO O POTENCIAL NO PONTO
            {
                double Haux = -ccc;
                ierr = MatSetValues(C, 1, &isp, 1, &isp, &Haux, ADD_VALUES);
            }
        }
    }

    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    //  MatView(C, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    // double val;
    // for (int i = 0; i < collocPoints_.size(); i++)
    // {
    //     for (int j = 0; j < collocPoints_.size(); j++)
    //     {
    //         ierr = MatGetValues(A, 1, &i, 1, &j, &val);
    //         CHKERRQ(ierr);
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    ///FAZER LOOPING ADICIONANDO O VALOR DO QUE A GENTE CONHECE NO PONTO DE COLOCAÇÃO NO VETOR x (CONDIÇÕES DE CONTORNO)
    double value;
    if (rank == 0)
    {
        for (int ic = 0; ic < ndegree; ic++)
        {
            if (collocCondition_[ic] == 0) ///CONHEÇO O FLUXO NO PONTO
            {
                value = collocPoints_[ic]->getFlux();
                ierr = VecSetValues(x, 1, &ic, &value, ADD_VALUES);
            }
            else ///CONHEÇO O POTENCIAL NO PONTO
            {
                value = collocPoints_[ic]->getPotential();
                ierr = VecSetValues(x, 1, &ic, &value, ADD_VALUES);
            }
        }
    }

    //Assemble vectors
    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    //Criando vetor do sistema
    ierr = MatMult(C, x, b);
    CHKERRQ(ierr);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    ////Solve using GMRES
    ierr = KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, 500);
    CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, PCBJACOBI);
    ierr = KSPSetType(ksp, KSPGMRES);
    CHKERRQ(ierr);
    ierr = KSPGMRESSetRestart(ksp, 500);
    CHKERRQ(ierr);

    //Solve using MUMPS
    // #if defined(PETSC_HAVE_MUMPS)
    //     ierr = KSPSetType(ksp, KSPPREONLY);
    //     ierr = KSPGetPC(ksp, &pc);
    //     ierr = PCSetType(pc, PCLU);
    // #endif

    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve linear system
    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(x, &ctx, &All);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    CHKERRQ(ierr);

    for (int ic = 0; ic < ndegree; ic++)
    {
        if (collocCondition_[ic] == 0) ///CONHEÇO O FLUXO NO PONTO
        {
            ierr = VecGetValues(All, 1, &ic, &value);
            CHKERRQ(ierr);
            collocPoints_[ic]->setPotential(value);
        }
        else ///CONHEÇO O POTENCIAL NO PONTO
        {
            ierr = VecGetValues(All, 1, &ic, &value);
            CHKERRQ(ierr);
            collocPoints_[ic]->setFlux(value);
        }
    }

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    if (internalPoints_.size() > 0)
    {
        computeInternalPointsPotentialProblem();
    }
    if (rank == 0)
    {
        exportToParaviewGeometricMesh(1);
        exportToParaviewSourcePoints(1);
        exportToParaviewCollocationMesh_Potential(1);
    }
    return 0;
}

void Problem::applyPotentialBoundaryConditions()
{
    collocCondition_.resize(collocPoints_.size(), 0); //todas as condições iniciais são de fluxo (Neumann)
    std::unordered_map<std::string, vecDouble> cond;
    std::string name;
    vecDouble values;

    cond = geometry_->getDirichletCondition();
    for (auto &cond : cond)
    {
        name = cond.first;
        values = cond.second;
        if (name[0] == 'l')
        {
            for (BoundaryElement *el : lineElements_[name])
            {
                for (CollocationPoint *coloc : el->getCollocationConnection())
                {
                    coloc->setPotential(values[0]);
                    collocCondition_[coloc->getIndex()] = 1;
                }
            }
        }
    }

    cond = geometry_->getNeumannCondition();
    for (auto &cond : cond)
    {
        name = cond.first;
        values = cond.second;
        if (name[0] == 'l')
        {
            for (BoundaryElement *el : lineElements_[name])
            {
                for (CollocationPoint *coloc : el->getCollocationConnection())
                {
                    coloc->setFlux(values[0]);
                }
            }
        }
    }
}

void Problem::addInternalPoints(std::vector<std::vector<double>> coord)
{
    for (int i = 0, iN = coord.size(); i < iN; i++)
    {
        SourcePoint *source = new SourcePoint(i, coord[i]);
        internalPoints_.push_back(source);
        internalPotential_.push_back(0.0);
        internalFlux_.push_back({0.0, 0.0});
    }
}

void Problem::addInternalPoints(Surface *surface, const std::vector<std::vector<double>> &coord)
{
    std::string surfaceName = surface->getName();
    if (subInternalPoints_.count(surfaceName) == 0)
    {
        std::vector<SourcePoint *> points;
        subInternalPoints_[surfaceName] = points;
    }

    int cont = internalPoints_.size();
    for (auto &points : coord)
    {
        SourcePoint *source = new SourcePoint(cont++, points);
        internalPoints_.push_back(source);
        subInternalPoints_[surfaceName].push_back(source);
        internalDisplacements_.push_back({0.0, 0.0});
        internalStresses_.push_back({0.0, 0.0, 0.0, 0.0});
    }
}

int Problem::computeInternalPointsPotentialProblem()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int ndegree = collocPoints_.size();
    const int nInterP = internalPoints_.size();

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          nInterP, ndegree, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, ndegree);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &All);
    CHKERRQ(ierr);
    ierr = VecSetSizes(All, PETSC_DECIDE, nInterP);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(All);
    CHKERRQ(ierr);

    matDouble matGl, matHl;
    if (rank == 0)
    {
        for (BoundaryElement *el : elements_)
        {
            el->potentialContribution(internalPoints_, matGl, matHl);
            std::vector<CollocationPoint *> conec = el->getCollocationConnection();

            for (int in = 0, nNodes = conec.size(); in < nNodes; in++)
            {
                int indexCP = conec[in]->getIndex();
                for (int isp = 0; isp < nInterP; isp++)
                {
                    ierr = MatSetValues(A, 1, &isp, 1, &indexCP, &matHl[isp][in], ADD_VALUES);
                    ierr = MatSetValues(C, 1, &isp, 1, &indexCP, &matGl[isp][in], ADD_VALUES);
                }
            }
        }
    }

    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    double value;
    if (rank == 0)
    {
        for (int ic = 0; ic < ndegree; ic++)
        {
            value = collocPoints_[ic]->getFlux();
            ierr = VecSetValues(x, 1, &ic, &value, ADD_VALUES);

            value = collocPoints_[ic]->getPotential();
            ierr = VecSetValues(b, 1, &ic, &value, ADD_VALUES);
        }
    }

    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);

    ierr = MatMult(C, x, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterP; ic++)
    {
        ierr = VecGetValues(All, 1, &ic, &value);
        CHKERRQ(ierr);
        internalPotential_[ic] = value;
    }

    ierr = MatMult(A, b, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterP; ic++)
    {
        ierr = VecGetValues(All, 1, &ic, &value);
        CHKERRQ(ierr);
        internalPotential_[ic] -= value;
    }

    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    //////////
    // FLUX //
    //////////

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          2 * nInterP, ndegree, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &All);
    CHKERRQ(ierr);
    ierr = VecSetSizes(All, PETSC_DECIDE, 2 * nInterP);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(All);
    CHKERRQ(ierr);

    if (rank == 0)
    {
        for (BoundaryElement *el : elements_)
        {
            el->calculateInternalFlux(internalPoints_, matGl, matHl);
            std::vector<CollocationPoint *> conec = el->getCollocationConnection();

            for (int in = 0, nNodes = conec.size(); in < nNodes; in++)
            {
                int indexCP = conec[in]->getIndex();
                for (int isp = 0; isp < 2 * nInterP; isp++)
                {
                    ierr = MatSetValues(A, 1, &isp, 1, &indexCP, &matHl[isp][in], ADD_VALUES);
                    ierr = MatSetValues(C, 1, &isp, 1, &indexCP, &matGl[isp][in], ADD_VALUES);
                }
            }
        }
    }
    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);
    // MatView(C, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    double val1, val2;
    ierr = MatMult(C, x, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterP; ic++)
    {
        Idof = 2 * ic;
        ierr = VecGetValues(All, 1, &Idof, &val1);
        CHKERRQ(ierr);

        Idof = 2 * ic + 1;
        ierr = VecGetValues(All, 1, &Idof, &val2);
        CHKERRQ(ierr);

        internalFlux_[ic] = {val1, val2};
    }

    ierr = MatMult(A, b, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterP; ic++)
    {
        Idof = 2 * ic;
        ierr = VecGetValues(All, 1, &Idof, &val1);
        CHKERRQ(ierr);

        Idof = 2 * ic + 1;
        ierr = VecGetValues(All, 1, &Idof, &val2);
        CHKERRQ(ierr);

        internalFlux_[ic][0] -= val1;
        internalFlux_[ic][1] -= val2;
        std::cout << internalPotential_[ic] << " " << internalFlux_[ic][0] << " " << internalFlux_[ic][1] << std::endl;
    }

    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    return 0;
}

int Problem::computeInternalPointsElasticityProblem()
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    const int nCollocationPoints = collocPoints_.size();
    const int nInterPoints = internalPoints_.size();

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          2 * nInterPoints, 2 * nCollocationPoints, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, 2 * nCollocationPoints);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &All);
    CHKERRQ(ierr);
    ierr = VecSetSizes(All, PETSC_DECIDE, 2 * nInterPoints);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(All);
    CHKERRQ(ierr);

    if (rank == 0)
    {
        matDouble matGl, matHl;
        for (auto &subRegionIntPoints : subInternalPoints_)
        {
            const int indexMat = geometry_->getIndexMaterial(subRegionIntPoints.first);
            const int nIntPointSubRegion = subRegionIntPoints.second.size();

            int internalPointIndex[2 * nIntPointSubRegion];
            for (int is = 0; is < nIntPointSubRegion; is++)
            {
                internalPointIndex[2 * is] = 2 * subRegionIntPoints.second[is]->getIndex();
                internalPointIndex[2 * is + 1] = 2 * subRegionIntPoints.second[is]->getIndex() + 1;
            }

            for (auto &el : subElements_[subRegionIntPoints.first])
            {
                el->elasticityContribution(subRegionIntPoints.second, matGl, matHl, materials_[indexMat]);
                std::vector<CollocationPoint *> conec = el->getCollocationConnection();

                for (int in = 0, nNodes = conec.size(); in < nNodes; in++)
                {
                    for (int dir = 0; dir < 2; dir++)
                    {
                        const int index = 2 * conec[in]->getIndex() + dir;
                        const int indexLocal = 2 * in + dir;
                        for (int iIP = 0; iIP < 2 * nIntPointSubRegion; iIP++)
                        {
                            ierr = MatSetValues(A, 1, &internalPointIndex[iIP], 1, &index, &matHl[iIP][indexLocal], ADD_VALUES);
                            ierr = MatSetValues(C, 1, &internalPointIndex[iIP], 1, &index, &matGl[iIP][indexLocal], ADD_VALUES);
                        }
                    }
                }
            }
        }
    }

    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    double value;

    if (rank == 0)
    {
        for (int ic = 0; ic < nCollocationPoints; ic++)
        {
            for (int dir = 0; dir < 2; dir++)
            {
                const int index = 2 * collocPoints_[ic]->getIndex() + dir;

                value = collocPoints_[ic]->getForce(dir);
                ierr = VecSetValues(x, 1, &index, &value, ADD_VALUES);

                value = collocPoints_[ic]->getDisplacement(dir);
                ierr = VecSetValues(b, 1, &index, &value, ADD_VALUES);
            }
        }
    }

    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(b);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(b);
    CHKERRQ(ierr);

    ierr = MatMult(C, x, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterPoints; ic++)
    {
        for (int dir = 0; dir < 2; dir++)
        {
            const int index = 2 * internalPoints_[ic]->getIndex() + dir;
            ierr = VecGetValues(All, 1, &index, &value);
            CHKERRQ(ierr);
            internalDisplacements_[ic][dir] = value;
        }
    }

    ierr = MatMult(A, b, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterPoints; ic++)
    {
        for (int dir = 0; dir < 2; dir++)
        {
            const int index = 2 * internalPoints_[ic]->getIndex() + dir;
            ierr = VecGetValues(All, 1, &index, &value);
            CHKERRQ(ierr);
            internalDisplacements_[ic][dir] -= value;
        }
    }

    if (bodyForces_.size() > 0)
    {
        vecDouble vecB;
        for (auto &bf : bodyForces_)
        {
            const std::string surfaceName = bf.first;
            const int indexMat = geometry_->getIndexMaterial(surfaceName);
            const std::vector<SourcePoint *> intPoints = subInternalPoints_[surfaceName];

            for (auto &el : subElements_[surfaceName])
            {
                el->elasticityBodyForceContribution(intPoints, vecB, materials_[indexMat], bf.second);
                for (int ip = 0, nIP = intPoints.size(); ip < nIP; ip++)
                {
                    int indexIP = intPoints[ip]->getIndex();
                    internalDisplacements_[indexIP][0] += vecB[2 * ip];
                    internalDisplacements_[indexIP][1] += vecB[2 * ip + 1];
                }
            }
        }
    }

    cout << "\nINTERNAL DISPLACEMENTS (X, Y):\n";
    for (int i = 0; i < nInterPoints; i++)
    {
        cout << internalDisplacements_[i][0] << " " << internalDisplacements_[i][1] << "\n";
    }

    ///////////////////
    // CAUCHY STRESS //
    ///////////////////

    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          4 * nInterPoints, 2 * nCollocationPoints, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);
    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &All);
    CHKERRQ(ierr);
    ierr = VecSetSizes(All, PETSC_DECIDE, 4 * nInterPoints);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(All);
    CHKERRQ(ierr);

    if (rank == 0)
    {
        matDouble matGl, matHl;
        for (auto &subRegionIntPoints : subInternalPoints_)
        {
            const int indexMat = geometry_->getIndexMaterial(subRegionIntPoints.first);
            const int nIntPointSubRegion = subRegionIntPoints.second.size();

            int internalPointIndex[4 * nIntPointSubRegion];
            for (int is = 0; is < nIntPointSubRegion; is++)
            {
                internalPointIndex[4 * is] = 4 * subRegionIntPoints.second[is]->getIndex();
                internalPointIndex[4 * is + 1] = 4 * subRegionIntPoints.second[is]->getIndex() + 1;
                internalPointIndex[4 * is + 2] = 4 * subRegionIntPoints.second[is]->getIndex() + 2;
                internalPointIndex[4 * is + 3] = 4 * subRegionIntPoints.second[is]->getIndex() + 3;
            }

            for (auto &el : subElements_[subRegionIntPoints.first])
            {
                el->calculateInternalStress(subRegionIntPoints.second, matGl, matHl, materials_[indexMat]);
                std::vector<CollocationPoint *> conec = el->getCollocationConnection();

                for (int in = 0, nNodes = conec.size(); in < nNodes; in++)
                {
                    for (int dir = 0; dir < 2; dir++)
                    {
                        const int index = 2 * conec[in]->getIndex() + dir;
                        const int indexLocal = 2 * in + dir;
                        for (int iIP = 0; iIP < 4 * nIntPointSubRegion; iIP++)
                        {
                            ierr = MatSetValues(A, 1, &internalPointIndex[iIP], 1, &index, &matHl[iIP][indexLocal], ADD_VALUES);
                            ierr = MatSetValues(C, 1, &internalPointIndex[iIP], 1, &index, &matGl[iIP][indexLocal], ADD_VALUES);
                        }
                    }
                }
            }
        }
    }

    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);
    // MatView(C, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    ierr = MatMult(C, x, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterPoints; ic++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            const int index = 4 * internalPoints_[ic]->getIndex() + dir;
            ierr = VecGetValues(All, 1, &index, &value);
            CHKERRQ(ierr);
            internalStresses_[ic][dir] = value;
        }
    }

    // VecView(All, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    ierr = MatMult(A, b, All);
    CHKERRQ(ierr);
    for (int ic = 0; ic < nInterPoints; ic++)
    {
        for (int dir = 0; dir < 4; dir++)
        {
            const int index = 4 * internalPoints_[ic]->getIndex() + dir;
            ierr = VecGetValues(All, 1, &index, &value);
            CHKERRQ(ierr);
            internalStresses_[ic][dir] -= value;
        }
    }

    if (bodyForces_.size() > 0)
    {
        vecDouble vecB;
        for (auto &bf : bodyForces_)
        {
            const std::string surfaceName = bf.first;
            const int indexMat = geometry_->getIndexMaterial(surfaceName);
            const std::vector<SourcePoint *> intPoints = subInternalPoints_[surfaceName];

            for (auto &el : subElements_[surfaceName])
            {
                el->calculateInternalStressBodyForce(intPoints, vecB, materials_[indexMat], bf.second);
                for (int ip = 0, nIP = intPoints.size(); ip < nIP; ip++)
                {
                    int indexIP = intPoints[ip]->getIndex();
                    internalStresses_[indexIP][0] += vecB[4 * ip];
                    internalStresses_[indexIP][1] += vecB[4 * ip + 1];
                    internalStresses_[indexIP][2] += vecB[4 * ip + 2];
                    internalStresses_[indexIP][3] += vecB[4 * ip + 3];
                }
            }
        }
    }

    // VecView(All, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    cout << "\nINTERNAL STRESS (XX, XY, YX, YY):\n";
    for (int i = 0; i < nInterPoints; i++)
    {
        cout << internalStresses_[i][0] << " " << internalStresses_[i][1] << " " << internalStresses_[i][2] << " " << internalStresses_[i][3] << "\n";
    }
    return 0;
}

void Problem::teste()
{
    // SourcePoint *s = sourcePoints_[16];
    // std::cout << s->getCoordinate(0) << " " << s->getCoordinate(1) << std::endl;

    // for (int i = 0; i < elements_.size(); i++)
    // {
    //     if (s->checkElement(i))
    //     {
    //         std::cout << "ELEMENTO " << i << " CONTEM...\n";
    //     }
    // }
    for (int i = 0; i < 4; i++)
    {

        std::string name = "l" + std::to_string(i);
        std::ofstream file(name + ".txt");
        for (BoundaryElement *el : lineElements_[name])
        {
            for (CollocationPoint *c : el->getCollocationConnection())
            {
                if (i == 0)
                {
                    file << c->getCoordinate(0) << " " << c->getPotential() << " " << c->getFlux() << "\n";
                }
                else if (i == 1)
                {
                    file << 4.0 + c->getCoordinate(1) << " " << c->getPotential() << " " << c->getFlux() << "\n";
                }
                else if (i == 2)
                {
                    file << 8.0 + (4.0 - c->getCoordinate(0)) << " " << c->getPotential() << " " << c->getFlux() << "\n";
                }
                else
                {
                    file << 12.0 + (4.0 - c->getCoordinate(1)) << " " << c->getPotential() << " " << c->getFlux() << "\n";
                }
            }
        }
    }
}

int Problem::solveElasticityProblem(const std::string &planeState)
{
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    applyElasticityBoundaryConditions();
    createSourcePoints();

    if (rank == 0)
    {
        exportToParaviewGeometricMesh(0);
        exportToParaviewSourcePoints(0);
        exportToParaviewCollocationMesh_Elasticity(0);
    }

    MPI_Barrier(PETSC_COMM_WORLD);

    int nCollocation = collocPoints_.size();
    int nDegree = 2 * nCollocation;

    //Create PETSc dense parallel matrix
    ierr = MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                          nDegree, nDegree, NULL, &A);
    CHKERRQ(ierr);
    ierr = MatGetOwnershipRange(A, &Istart, &Iend);
    CHKERRQ(ierr);
    ierr = MatDuplicate(A, MAT_DO_NOT_COPY_VALUES, &C);
    CHKERRQ(ierr);

    //Create PETSc vectors
    ierr = VecCreate(PETSC_COMM_WORLD, &b);
    CHKERRQ(ierr);
    ierr = VecSetSizes(b, PETSC_DECIDE, nDegree);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(b);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &x);
    CHKERRQ(ierr);
    ierr = VecDuplicate(b, &All);
    CHKERRQ(ierr);

    matDouble matGl, matHl;
    if (rank == 0)
    {
        for (auto &sub : subElements_)
        {
            const int indexMat = geometry_->getIndexMaterial(sub.first);
            const std::vector<SourcePoint *> sourcePoints = subSourcePoints_[sub.first];
            const int nSourceP = sourcePoints.size();
            int indexSourcePoints[2 * nSourceP];
            for (int is = 0; is < nSourceP; is++)
            {
                indexSourcePoints[2 * is] = 2 * sourcePoints[is]->getIndex();
                indexSourcePoints[2 * is + 1] = 2 * sourcePoints[is]->getIndex() + 1;
            }

            for (auto &el : sub.second)
            {
                el->elasticityContribution(sourcePoints, matGl, matHl, materials_[indexMat]);
                std::vector<CollocationPoint *> conec = el->getCollocationConnection();
                for (int ic = 0, nConec = conec.size(); ic < nConec; ic++)
                {
                    for (int dir = 0; dir < 2; dir++)
                    {
                        const int index = 2 * conec[ic]->getIndex() + dir;
                        const int indexLocal = 2 * ic + dir;
                        const int cond = collocCondition_[index];
                        if (cond == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
                        {
                            for (int isp = 0; isp < 2 * nSourceP; isp++)
                            {
                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &index, &matHl[isp][indexLocal], ADD_VALUES); //A é a matriz H e C é a matriz G
                                ierr = MatSetValues(C, 1, &indexSourcePoints[isp], 1, &index, &matGl[isp][indexLocal], ADD_VALUES);
                            }
                        }
                        else if (cond == 1) ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
                        {
                            for (int isp = 0; isp < 2 * nSourceP; isp++)
                            {
                                double Haux = -matHl[isp][indexLocal];
                                double Gaux = -matGl[isp][indexLocal];

                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &index, &Gaux, ADD_VALUES);
                                ierr = MatSetValues(C, 1, &indexSourcePoints[isp], 1, &index, &Haux, ADD_VALUES);
                            }
                        }
                        else if (cond == 2)
                        {
                            const int indexAux = 2 * coupledCollocFirst_[conec[ic]->getIndex()] + dir;
                            for (int isp = 0; isp < 2 * nSourceP; isp++)
                            {
                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &index, &matHl[isp][indexLocal], ADD_VALUES); //A é a matriz H e C é a matriz G

                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &indexAux, &matGl[isp][indexLocal], ADD_VALUES);
                            }
                        }
                        else if (cond == 3)
                        {
                            const int indexAux = 2 * coupledCollocSecond_[conec[ic]->getIndex()] + dir;
                            for (int isp = 0; isp < 2 * nSourceP; isp++)
                            {
                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &indexAux, &matHl[isp][indexLocal], ADD_VALUES); //A é a matriz H e C é a matriz G

                                double Gaux = -matGl[isp][indexLocal];
                                ierr = MatSetValues(A, 1, &indexSourcePoints[isp], 1, &index, &Gaux, ADD_VALUES);
                            }
                        }
                        else
                        {
                            cout << "ERROR\n";
                        }
                    }
                }
            }
        }
    }

    if (!sourceOut_ and rank == 0)
    {
        double ccc = 0.5;
        for (int icoloc = 0; icoloc < nCollocation; icoloc++)
        {
            for (int dir = 0; dir < 2; dir++)
            {
                const int isp = 2 * icoloc + dir;
                const int cond = collocCondition_[isp];
                if (cond == 0) ///CONHEÇO A FORÇA NO PONTO//DIREÇÃO
                {
                    ierr = MatSetValues(A, 1, &isp, 1, &isp, &ccc, ADD_VALUES);
                }
                else if (cond == 1) ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
                {
                    double Haux = -ccc;
                    ierr = MatSetValues(C, 1, &isp, 1, &isp, &Haux, ADD_VALUES);
                }
                else if (cond == 2)
                {
                    ierr = MatSetValues(A, 1, &isp, 1, &isp, &ccc, ADD_VALUES);
                }
                else if (cond == 3)
                {
                    const int indexAux = 2 * coupledCollocSecond_[icoloc] + dir;
                    ierr = MatSetValues(A, 1, &isp, 1, &indexAux, &ccc, ADD_VALUES);
                }
            }
        }
    }

    //Assemble matrices
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    //  MatView(C, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    // double val;
    // for (int i = 0; i < collocPoints_.size(); i++)
    // {
    //     for (int j = 0; j < collocPoints_.size(); j++)
    //     {
    //         ierr = MatGetValues(A, 1, &i, 1, &j, &val);
    //         CHKERRQ(ierr);
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    ///FAZER LOOPING ADICIONANDO O VALOR DO QUE A GENTE CONHECE NO PONTO DE COLOCAÇÃO NO VETOR x (CONDIÇÕES DE CONTORNO)
    double value;
    if (rank == 0)
    {
        for (int ic = 0; ic < nCollocation; ic++)
        {
            for (int dir = 0; dir < 2; dir++)
            {
                const int index = 2 * ic + dir;
                const int cond = collocCondition_[index];
                if (cond == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
                {
                    value = collocPoints_[ic]->getForce(dir);
                    ierr = VecSetValues(x, 1, &index, &value, ADD_VALUES);
                }
                else if (cond == 1) ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
                {
                    value = collocPoints_[ic]->getDisplacement(dir);
                    ierr = VecSetValues(x, 1, &index, &value, ADD_VALUES);
                }
                // else if (cond == 2)
                // {
                // }
                // else if (cond == 3)
                // {
                // }
            }
        }
    }

    //Assemble vectors
    ierr = VecAssemblyBegin(x);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    CHKERRQ(ierr);

    // VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    //Criando vetor do sistema
    ierr = MatMult(C, x, b);
    CHKERRQ(ierr);

    if (bodyForces_.size() > 0)
    {
        vecDouble vecB;
        for (auto &bf : bodyForces_)
        {
            std::string surfaceName = bf.first;
            const int indexMat = geometry_->getIndexMaterial(surfaceName);
            const std::vector<SourcePoint *> sourcePoints = subSourcePoints_[surfaceName];
            const int nSourceP = sourcePoints.size();
            int indexSourcePoints[2 * nSourceP];
            for (int is = 0; is < nSourceP; is++)
            {
                indexSourcePoints[2 * is] = 2 * sourcePoints[is]->getIndex();
                indexSourcePoints[2 * is + 1] = 2 * sourcePoints[is]->getIndex() + 1;
            }

            for (auto &el : subElements_[surfaceName])
            {
                el->elasticityBodyForceContribution(sourcePoints, vecB, materials_[indexMat], bf.second);
                for (int isp = 0; isp < 2 * nSourceP; isp++)
                {
                    ierr = VecSetValues(b, 1, &indexSourcePoints[isp], &vecB[isp], ADD_VALUES);
                }
            }
        }
        //Assemble vectors
        ierr = VecAssemblyBegin(b);
        CHKERRQ(ierr);
        ierr = VecAssemblyEnd(b);
        CHKERRQ(ierr);
    }

    // VecView(b, PETSC_VIEWER_STDOUT_WORLD);
    // CHKERRQ(ierr);

    //Create KSP context to solve the linear system
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A);
    CHKERRQ(ierr);

    //Solve using GMRES
    // ierr = KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
    // CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp, &pc);
    // ierr = PCSetType(pc, PCBJACOBI);
    // ierr = KSPSetType(ksp, KSPGMRES);
    // CHKERRQ(ierr);
    // ierr = KSPGMRESSetRestart(ksp, 1000);
    // CHKERRQ(ierr);

    // Solve using MUMPS
#if defined(PETSC_HAVE_MUMPS)
    ierr = KSPSetType(ksp, KSPPREONLY);
    ierr = KSPGetPC(ksp, &pc);
    ierr = PCSetType(pc, PCLU);
#endif

    ierr = KSPSetFromOptions(ksp);
    CHKERRQ(ierr);
    ierr = KSPSetUp(ksp);

    //Solve linear system
    ierr = KSPSolve(ksp, b, x);
    CHKERRQ(ierr);
    ierr = KSPGetTotalIterations(ksp, &iterations);

    //Gathers the solution vector to the master process
    ierr = VecScatterCreateToAll(x, &ctx, &All);
    CHKERRQ(ierr);
    ierr = VecScatterBegin(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterEnd(ctx, x, All, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
    ierr = VecScatterDestroy(&ctx);
    CHKERRQ(ierr);

    for (int ic = 0; ic < nCollocation; ic++)
    {
        for (int dir = 0; dir < 2; dir++)
        {
            const int index = 2 * ic + dir;
            const int cond = collocCondition_[index];
            if (cond == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setDisplacement(value, dir);
            }
            else if (cond == 1) ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setForce(value, dir);
            }
            else if (cond == 2)
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setDisplacement(value, dir);
                collocPoints_[coupledCollocFirst_[ic]]->setDisplacement(value, dir);
            }
            else if (cond == 3)
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setForce(value, dir);
                collocPoints_[coupledCollocSecond_[ic]]->setForce(-value, dir);
            }
        }
    }

    ierr = KSPDestroy(&ksp);
    CHKERRQ(ierr);
    ierr = VecDestroy(&b);
    CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    CHKERRQ(ierr);
    ierr = VecDestroy(&All);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&C);
    CHKERRQ(ierr);

    if (internalPoints_.size() > 0)
    {
        computeInternalPointsElasticityProblem();
    }
    if (rank == 0)
    {
        exportToParaviewGeometricMesh(1);
        exportToParaviewSourcePoints(1);
        exportToParaviewCollocationMesh_Elasticity(1);
    }
    return 0;
}

void Problem::applyElasticityBoundaryConditions()
{
    collocCondition_.resize(2 * collocPoints_.size(), 0); //todas as condições iniciais são de força (Neumann)
    std::unordered_map<std::string, vecDouble> cond;
    std::string name;
    vecDouble values;

    cond = geometry_->getDirichletCondition();
    for (auto &cond : cond)
    {
        name = cond.first;
        values = cond.second;
        if (name[0] == 'l')
        {
            for (BoundaryElement *el : lineElements_[name])
            {
                for (CollocationPoint *coloc : el->getCollocationConnection())
                {
                    int index = coloc->getIndex();
                    if (values[0] != 1.0e-240)
                    {
                        coloc->setDisplacement(values[0], 0);
                        collocCondition_[2 * index] = 1;
                    }
                    if (values[1] != 1.0e-240)
                    {
                        coloc->setDisplacement(values[1], 1);
                        collocCondition_[2 * index + 1] = 1;
                    }
                }
            }
        }
    }

    cond = geometry_->getNeumannCondition();
    for (auto &cond : cond)
    {
        name = cond.first;
        values = cond.second;
        if (name[0] == 'l')
        {
            for (BoundaryElement *el : lineElements_[name])
            {
                for (CollocationPoint *coloc : el->getCollocationConnection())
                {
                    int index = coloc->getIndex();
                    if (values[0] != 1.0e-240)
                    {
                        coloc->setForce(values[0], 0);
                    }
                    if (values[1] != 1.0e-240)
                    {
                        coloc->setForce(values[1], 1);
                    }
                }
            }
        }
    }

    for (auto &c : coupledCollocFirst_)
    {
        collocCondition_[2 * c.first] = 2;
        collocCondition_[2 * c.first + 1] = 2;

        collocCondition_[2 * c.second] = 3;
        collocCondition_[2 * c.second + 1] = 3;
    }
}

void Problem::exportToParaviewCollocationMesh_Elasticity(const int &index)
{
    std::stringstream nameFile;
    nameFile << "collocationMesh" << index << ".vtu";
    std::ofstream file(nameFile.str());

    int nElement = 0;
    for (auto &el : lineElements_)
    {
        nElement += el.second.size();
    }

    //header
    file << "<?xml version=\"1.0\"?>"
         << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">"
         << "\n"
         << "  <UnstructuredGrid>"
         << "\n"
         << "  <Piece NumberOfPoints=\"" << collocPoints_.size() + internalPoints_.size()
         << "\"  NumberOfCells=\"" << nElement + internalPoints_.size()
         << "\">"
         << "\n";
    //nodal coordinates
    file << "    <Points>"
         << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"3\" format=\"ascii\">"
         << "\n";

    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getCurrentCoordinate(0) << " " << n->getCurrentCoordinate(1) << " " << n->getCurrentCoordinate(2) << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << n->getCoordinate(0) + internalDisplacements_[n->getIndex()][0] << " " << n->getCoordinate(1) + internalDisplacements_[n->getIndex()][1] << " " << n->getCoordinate(2) << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Points>"
         << "\n";
    //element connectivity
    file << "    <Cells>"
         << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">"
         << "\n";

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            for (Node *no : el->getConnection())
            {
                file << no->getIndex() << " ";
            }
            file << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << collocPoints_.size() + n->getIndex() << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";
    int aux = 0;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            aux += el->getConnection().size();
            file << aux << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << ++aux << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 68 << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << 1 << "\n";
    }

    file << "      </DataArray>"
         << "\n"
         << "    </Cells>"
         << "\n";
    //nodal results
    file << "    <PointData>"
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Displacement\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getDisplacement(0) << " " << n->getDisplacement(1) << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << internalDisplacements_[n->getIndex()][0] << " " << internalDisplacements_[n->getIndex()][1] << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Traction\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getForce(0) << " " << n->getForce(1) << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << 0.0 << " " << 0.0 << "\n";
        // file << sqrt(internalFlux_[n->getIndex()][0] * internalFlux_[n->getIndex()][0] + internalFlux_[n->getIndex()][1] * internalFlux_[n->getIndex()][1]) << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    // file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
    //      << "Name=\"InternalFlux\" format=\"ascii\">"
    //      << "\n";
    // for (CollocationPoint *n : collocPoints_)
    // {
    //     file << 0 << " " << 0 << "\n";
    // }
    // for (SourcePoint *n : internalPoints_)
    // {
    //     file << internalFlux_[n->getIndex()][0] << " " << internalFlux_[n->getIndex()][1] << "\n";
    // }
    // file << "      </DataArray> "
    //      << "\n";

    file << "    </PointData>"
         << "\n";

    //elemental results
    file << "    <CellData>"
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"Process\" format=\"ascii\">" << std::endl;

    for (int ie = 0; ie < nElement; ie++)
    {
        file << 0 << "\n";
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << 0 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"LineNumber\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        std::string ln = line.first.substr(1, line.first.size() - 1);
        for (BoundaryElement *el : line.second)
        {
            file << ln << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << -1 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "<DataArray type=\"Float64\" NumberOfComponents=\"1\" "
         << "Name=\"ElementIndex\" format=\"ascii\">" << std::endl;

    for (auto &line : lineElements_)
    {
        for (BoundaryElement *el : line.second)
        {
            file << el->getIndex() << "\n";
        }
    }
    for (SourcePoint *n : internalPoints_)
    {
        file << -1 << "\n";
    }
    file << "      </DataArray> "
         << "\n";

    file << "    </CellData>"
         << "\n";
    //footnote
    file << "  </Piece>"
         << "\n"
         << "  </UnstructuredGrid>"
         << "\n"
         << "</VTKFile>"
         << "\n";
}

void Problem::addMaterial(const double &young, const double &poisson, const double &density)
{
    Material *mat = new Material(young, poisson, density);
    materials_.push_back(mat);
}

void Problem::coupleLines(std::vector<std::vector<Line *>> coupledLines)
{
    for (auto &c : coupledLines)
    {
        if (c.size() == 2)
        {
            std::unordered_set<CollocationPoint *> line0, line1;
            for (auto &el : lineElements_[c[0]->getName()])
            {
                for (auto &colloc : el->getCollocationConnection())
                {
                    line0.insert(colloc);
                }
            }
            for (auto &el : lineElements_[c[1]->getName()])
            {
                for (auto &colloc : el->getCollocationConnection())
                {
                    line1.insert(colloc);
                }
            }
            vecDouble coord0(3), coord1(3);
            for (auto &colloc0 : line0)
            {
                bool find = false;
                colloc0->getCoordinate(coord0);
                for (auto &colloc1 : line1)
                {
                    colloc1->getCoordinate(coord1);
                    if (fabs(coord0[0] - coord1[0]) <= 1.0e-06 and fabs(coord0[1] - coord1[1]) <= 1.0e-06)
                    {
                        coupledCollocFirst_[colloc0->getIndex()] = colloc1->getIndex();
                        coupledCollocSecond_[colloc1->getIndex()] = colloc0->getIndex();
                        find = true;
                        break;
                    }
                }
                if (!find)
                {
                    cout << "COLOCAÇÃO " << colloc0->getIndex() << " NÃO ENCONTROU SEU PAR!\n";
                }
            }
        }
        else
        {
            cout << "CHECAR LINHAS ACOPLADAS\n";
        }
    }
    // for (auto &tt : coupledCollocFirst_)
    // {
    //     cout << tt.first << " " << tt.second << "\n";
    // }
}

void Problem::addBodyForceToSurface(Surface *surface, const vecDouble &values)
{
    bodyForces_[surface->getName()] = values;
}

void Problem::teste2()
{
    SourcePoint *sp = new SourcePoint(0, {0.5, 0.5});

    vecDouble vecB;
    vecDouble result(2, 0.0);
    for (BoundaryElement *el : elements_)
    {
        el->elasticityBodyForceContribution({sp}, vecB, materials_[0], {1.0, 0.0});
        result[0] += vecB[0];
        result[1] += vecB[1];
    }
    cout << vecB.size() << ": " << result[0] << " " << result[1] << endl;
}
