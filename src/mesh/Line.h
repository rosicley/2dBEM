#pragma once

#include "Point.h"
#include "../Node.h"

class Line
{
public:
    Line();

    Line(const std::string &name, std::vector<Point *> points, const bool &discretization);

    virtual ~Line();

    Line *operator-();

    std::string getName();

    Point *getInitialPoint();

    Point *getEndPoint();

    virtual std::string getGmshCode();

    void setName(const std::string &name);

    void setPoints(const std::vector<Point *> &newPoints);

protected:
    std::string name_;
    std::vector<Point *> points_;
    bool discretization_;
};

class Circle : public Line
{
public:
    Circle(const std::string &name, std::vector<Point *> points, const bool &discretization);

    ~Circle();

    Circle *operator-();

    std::string getGmshCode() override;
};

class Spline : public Line
{
public:
    Spline(const std::string &name, std::vector<Point *> points, const bool &discretization = true);

    ~Spline();

    Spline *operator-();

    std::string getGmshCode() override;
};

class LineLoop
{
public:
    LineLoop();

    LineLoop(const std::string &name, std::vector<Line *> lines);

    ~LineLoop();

    LineLoop *operator-();

    std::string getName();

    Line *getLine(const int &index);

    std::vector<Line *> getLines();

    std::string getGmshCode();

    void verification();

private:
    std::string name_;
    std::vector<Line *> lines_;
};
