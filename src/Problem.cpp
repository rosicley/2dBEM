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

    /////CHECKING DUPLICATED NODES
    vecDouble coordRef(3), coordLoop(3);
    double tol = 1.0e-06;
    for (Node *noRef : nodes_)
    {
        noRef->getCoordinate(coordRef);
        for (Node *noLoop : nodes_)
        {
            noLoop->getCoordinate(coordLoop);
            if (fabs(coordRef[0] - coordLoop[0]) <= tol and fabs(coordRef[1] - coordLoop[1]) <= tol and noRef != noLoop)
            {
                duplicatedNodes_.insert(noRef);
            }
        }
    }
    //////CHECKING DUPLICATED NODES

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
    int node0 = duplicatedNodes_.count(connection[0]);
    int node1 = duplicatedNodes_.count(connection[1]);
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
            int node0 = duplicatedNodes_.count(connection[0]);
            int node1 = duplicatedNodes_.count(connection[1]);
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
        file << 4 << "\n";
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
         << "Name=\"IndexElement\" format=\"ascii\">" << std::endl;

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
        file << 4 << "\n";
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
         << "Name=\"IndexElement\" format=\"ascii\">" << std::endl;

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

    for (SourcePoint *n : sourcePoints_)
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
    int cont = 0;
    for (SourcePoint *n : sourcePoints_)
    {
        file << cont++ << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //offsets
    file << "      <DataArray type=\"Int32\""
         << " Name=\"offsets\" format=\"ascii\">"
         << "\n";

    cont = 0;
    for (SourcePoint *n : sourcePoints_)
    {
        file << cont++ << "\n";
    }

    file << "      </DataArray>"
         << "\n";
    //elements type
    file << "      <DataArray type=\"UInt8\" Name=\"types\" "
         << "format=\"ascii\">"
         << "\n";

    for (SourcePoint *n : sourcePoints_)
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
        computeInternalPoints();
    }
    if (rank == 0)
    {
        exportToParaviewGeometricMesh(1);
        exportToParaviewSourcePoints(1);
        exportToParaviewCollocationMesh_Potential(1);
    }
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

int Problem::computeInternalPoints()
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

    // for (int i = 0; i < nCollocation; i++)
    // {
    //     if (collocCondition_[2 * i] == 0)
    //     {
    //         cout << "FORCE: " << collocPoints_[i]->getForce(0);
    //     }
    //     else
    //     {
    //         cout << "DISPL: " << collocPoints_[i]->getDisplacement(0);
    //     }
    //     cout << " ";
    //     if (collocCondition_[2 * i + 1] == 0)
    //     {
    //         cout << "FORCE: " << collocPoints_[i]->getForce(1);
    //     }
    //     else
    //     {
    //         cout << "DISPL: " << collocPoints_[i]->getDisplacement(1);
    //     }
    //     cout << "\n";
    // }

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
        for (BoundaryElement *el : elements_)
        {
            el->elasticityContribution(sourcePoints_, matGl, matHl, materials_[0]);
            std::vector<CollocationPoint *> conec = el->getCollocationConnection();

            for (int ic = 0, nConec = conec.size(); ic < nConec; ic++)
            {
                for (int dir = 0; dir < 2; dir++)
                {
                    int index = 2 * conec[ic]->getIndex() + dir;
                    int indexLocal = 2 * ic + dir;
                    if (collocCondition_[index] == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
                    {
                        for (int isp = 0; isp < 2 * sourcePoints_.size(); isp++)
                        {
                            ierr = MatSetValues(A, 1, &isp, 1, &index, &matHl[isp][indexLocal], ADD_VALUES); //A é a matriz H e C é a matriz G
                            ierr = MatSetValues(C, 1, &isp, 1, &index, &matGl[isp][indexLocal], ADD_VALUES);
                        }
                    }
                    else ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
                    {
                        for (int isp = 0; isp < 2 * sourcePoints_.size(); isp++)
                        {
                            double Haux = -matHl[isp][indexLocal];
                            double Gaux = -matGl[isp][indexLocal];

                            ierr = MatSetValues(A, 1, &isp, 1, &index, &Gaux, ADD_VALUES);
                            ierr = MatSetValues(C, 1, &isp, 1, &index, &Haux, ADD_VALUES);
                        }
                    }
                }
            }
        }
    }

    if (!sourceOut_ and rank == 0)
    {
        double ccc = 0.5;
        for (int isp = 0; isp < nDegree; isp++)
        {
            if (collocCondition_[isp] == 0) ///CONHEÇO A FORÇA NO PONTO//DIREÇÃO
            {
                ierr = MatSetValues(A, 1, &isp, 1, &isp, &ccc, ADD_VALUES);
            }
            else ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
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
        for (int ic = 0; ic < nCollocation; ic++)
        {
            for (int dir = 0; dir < 2; dir++)
            {
                int index = 2 * ic + dir;
                if (collocCondition_[2 * ic + dir] == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
                {
                    value = collocPoints_[ic]->getForce(dir);
                    ierr = VecSetValues(x, 1, &index, &value, ADD_VALUES);
                }
                else ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
                {
                    value = collocPoints_[ic]->getDisplacement(dir);
                    ierr = VecSetValues(x, 1, &index, &value, ADD_VALUES);
                }
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
            int index = 2 * ic + dir;
            if (collocCondition_[2 * ic + dir] == 0) ///CONHEÇO A FORÇA NO PONTO/DIREÇÃO
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setDisplacement(value, dir);
            }
            else ///CONHEÇO O DESLOCAMENTO NO PONTO/DIREÇÃO
            {
                ierr = VecGetValues(All, 1, &index, &value);
                CHKERRQ(ierr);
                collocPoints_[ic]->setForce(value, dir);
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
        computeInternalPoints();
    }
    if (rank == 0)
    {
        exportToParaviewGeometricMesh(1);
        exportToParaviewSourcePoints(1);
        exportToParaviewCollocationMesh_Elasticity(1);
    }
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
    // for (SourcePoint *n : internalPoints_)
    // {
    //     file << n->getCoordinate(0) << " " << n->getCoordinate(1) << " " << n->getCoordinate(2) << "\n";
    // }

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
        file << 4 << "\n";
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
    // for (SourcePoint *n : internalPoints_)
    // {
    //     file << internalPotential_[n->getIndex()] << "\n";
    // }
    file << "      </DataArray> "
         << "\n";

    file << "      <DataArray type=\"Float64\" NumberOfComponents=\"2\" "
         << "Name=\"Traction\" format=\"ascii\">"
         << "\n";
    for (CollocationPoint *n : collocPoints_)
    {
        file << n->getForce(0) << " " << n->getForce(1) << "\n";
    }
    // for (SourcePoint *n : internalPoints_)
    // {
    //     file << sqrt(internalFlux_[n->getIndex()][0] * internalFlux_[n->getIndex()][0] + internalFlux_[n->getIndex()][1] * internalFlux_[n->getIndex()][1]) << "\n";
    // }
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
         << "Name=\"IndexElement\" format=\"ascii\">" << std::endl;

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
