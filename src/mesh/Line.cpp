#include "Line.h"

//////////
// Line //
//////////

Line::Line() {}

Line::Line(const std::string &name, std::vector<Point *> points, const bool &discretization)
{
    name_ = name;
    discretization_ = discretization;
    if (points.size() != 2)
    {
        cout << "A line must be defined with two points.\n";
        exit(EXIT_FAILURE);
    }
    points_.resize(2);
    points_[0] = points[0];
    points_[1] = points[1];
}

Line::~Line() {}

Line *Line::operator-()
{
    Line *copy = new Line("-" + name_, points_, discretization_);
    return copy;
}

std::string Line::getName()
{
    return name_;
}
Point *Line::getInitialPoint()
{
    return points_[0];
}

Point *Line::getEndPoint()
{
    return points_[points_.size() - 1];
}

std::string Line::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName();
    if (discretization_ == true)
    {

        text << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    }
    else
    {
        text << "};\n//\n";
    }
    return text.str();
}

void Line::setName(const std::string &name)
{
    name_ = name;
}

void Line::setPoints(const std::vector<Point *> &newPoints)
{
    points_ = newPoints;
}

////////////
// Circle //
////////////

Circle::Circle(const std::string &name, std::vector<Point *> points, const bool &discretization)
{
    name_ = name;
    discretization_ = discretization;
    if (points.size() != 3)
    {
        cout << "A circle must be defined with three points.\n";
        exit(EXIT_FAILURE);
    }
    points_.resize(3);
    points_[0] = points[0];
    points_[1] = points[1];
    points_[2] = points[2];
}

Circle::~Circle() {}

Circle *Circle::operator-()
{
    Circle *copy = new Circle(name_, points_, discretization_);
    copy->setName("-" + name_);
    return copy;
}

std::string Circle::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newl; Circle(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName() << ", " << points_[2]->getName();
    if (discretization_)
    {
        text << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    }
    else
    {
        text << "};\n//\n";
    }
    return text.str();
}

////////////
// Spline //
////////////

Spline::Spline(const std::string &name, std::vector<Point *> points, const bool &discretization)
{
    name_ = name;
    discretization_ = discretization;
    int nP = points.size();
    points_.resize(nP);
    for (int i = 0; i < nP; i++)
    {
        points_[i] = points[i];
    }
}

Spline::~Spline() {}

Spline *Spline::operator-()
{
    Spline *copy = new Spline(name_, points_, discretization_);
    copy->setName("-" + name_);
    return copy;
}

std::string Spline::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newl; Spline(" << name_ << ") = {";

    for (int i = 0, nPoints = points_.size(); i < nPoints; i++)
    {
        text << points_[i]->getName();
        if (i != nPoints - 1)
        {
            text << ", ";
        }
    }
    if (discretization_)
    {
        text << "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
    }
    else
    {
        text << "};\n//\n";
    }
    return text.str();
}

//////////////
// LineLoop //
//////////////

LineLoop::LineLoop() {}

LineLoop::LineLoop(const std::string &name, std::vector<Line *> lines)
{
    name_ = name;
    lines_.reserve(lines.size());
    for (Line *line : lines)
    {
        lines_.push_back(line);
    }
}

LineLoop::~LineLoop() {}

std::string LineLoop::getName()
{
    return name_;
}

LineLoop *LineLoop::operator-()
{
    LineLoop *copy = new LineLoop("-" + name_, lines_);
    return copy;
}

Line *LineLoop::getLine(const int &index)
{
    return lines_[index];
}

std::vector<Line *> LineLoop::getLines()
{
    return lines_;
}

std::string LineLoop::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newll; Line Loop(" << name_ << ") = {";
    for (int i = 0, nLine = lines_.size(); i < nLine; i++)
    {
        text << lines_[i]->getName();
        if (i != (nLine - 1))
            text << ", ";
    }
    text << "};\n//\n";
    return text.str();
}

void LineLoop::verification()
{
    for (int i = 0, nLine = lines_.size(); i < nLine; i++)
    {
        std::string name = lines_[i]->getName();
        int index = i;
        std::string previous_name = (i == 0) ? lines_[nLine - 1]->getName() : lines_[i - 1]->getName();
        int previous_index = (i == 0) ? nLine - 1 : i - 1;
        Point *initial_point = (name[0] == '-') ? lines_[index]->getEndPoint() : lines_[index]->getInitialPoint();
        Point *end_point = (previous_name[0] == '-') ? lines_[previous_index]->getInitialPoint() : lines_[previous_index]->getEndPoint();
        if (initial_point->getName() != end_point->getName())
        {
            std::stringstream text;
            text << "The lines ";
            for (int j = 0; j < nLine; j++)
            {
                text << lines_[j]->getName();
                if (j != (nLine - 1))
                    text << ", ";
            }
            text << " do not form a closed loop." << std::endl;
            std::cout << text.str();
            break;
        }
    }
}