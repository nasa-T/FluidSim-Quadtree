#include <cmath>
#include <iostream>
#include <fstream>
#include "fluid.h"

#include <map>
#include <string>
#include <functional>

#include <omp.h>
#include <algorithm>

#include <chrono>
#include <cassert>

class VelocityVector {
    public:
        VelocityVector(double vx_ = 0, double vy_ = 0, double ax_ = 0, double ay_ = 0): vx(vx_), vy(vy_), ax(ax_), ay(ay_) {}
        VelocityVector operator+(VelocityVector& other) {
            VelocityVector added(this->vx + other.vx, this->vy + other.vy, other.ax, other.ay);
            return added;
        }
        VelocityVector operator/(double d) {
            VelocityVector divided(this->vx/d, this->vy/d, ax, ay);
            return divided;
        }
        VelocityVector operator-(VelocityVector& other) {
            VelocityVector subtracted(this->vx - other.vx, this->vy - other.vy, other.ax, other.ay);
            return subtracted;
        }
        VelocityVector operator*(double d) {
            VelocityVector divided(this->vx*d, this->vy*d, ax, ay);
            return divided;
        }
        double getVx() {
            return vx;
        }
        double getVy() {
            return vy;
        }
        double getAx() {
            return ax;
        }
        double getAy() {
            return ay;
        }
        void setVx(double v) {
            vx = v;
        }
        void setVy(double v) {
            vy = v;
        }
        void accelerate() {
            vx += ax * consts::dt;
            vy += ay * consts::dt;
        }
        double getMag() {
            return std::sqrt(getVx()*getVx()+getVy()*getVy());
            // return pow(vx*vx+vy*vy,1/2);
        }
        // VelocityVector operator+(VelocityVector &vec) {
        //     return VelocityVector(vx+vec.getVx(), vy+vec.getVy());
        // }
        // std::ostream& operator<<(std::ostream &s, VelocityVector &vec) {
        //     return s << "vx: " << vec.getVx() << " vy: " << vec.getVy();
        // }
    private:
        double vx, vy, ax, ay;
};

struct Neighbors {
    FluidCell *top, *bottom, *left, *right = nullptr;
};

class FluidCell {
    public:
        FluidCell(FluidCell *parent, uint32_t depth, uint32_t idx, long double mass, double width, double height, double temp, float X=0.8, float Y=0.15, float Z=0.05): parent(parent), depth(depth), index(idx), mass(mass), width(width), height(height), temperature(temp), nw(nullptr), ne(nullptr), sw(nullptr), se(nullptr) {
            // size is the physical area
            size = width * height;
            density = mass/size;
            setX(X);
            setY(Y);
            setZ(Z);
            long double nonMet = (hydrogen+helium);
            long double cv = consts::r/(consts::HMol*hydrogen/nonMet+consts::HeMol*helium/nonMet)/(1.4-1);
            e = cv*temperature*density;
            // std::cout << e << std::endl;
            pressure = (1.4-1)*e;
            velocity = VelocityVector();
            neighbors = Neighbors();
            double gravPotential = 0;
            nw, ne, sw, se = nullptr;
        }

        long double getCv() {
            double nonMet = (hydrogen+helium);
            return consts::r/(consts::HMol*hydrogen/nonMet+consts::HeMol*helium/nonMet)/(1.4-1);
        }

        long double getMass() {
            return mass;
        }

        long double getSize() {
            return size;
        }

        long double getWidth() {
            return width;
        }

        long double getHeight() {
            return height;
        }

        long double getDensity() {
            return mass/size;
        }

        long double getTemp(bool recalculate=true) {
            if (recalculate) {
                // double cv = consts::r/consts::HMol/(1.4-1);
                if (mass == 0) {
                    // std::cout << mass << std::endl;
                    temperature = 0;
                } else {
                    temperature = gete(false)/(getDensity()*getCv());
                }
            }
            return temperature;
        }

        void setTemp(long double t) {
            temperature = t;
        }

        double getPressure(bool recalculate=true) {
            // double oldPressure = pressure;
            if (recalculate) {
                // pressure = getDensity()/consts::HMol * consts::r * getTemp();
                pressure = (1.4 - 1)*gete(false);
            }
            
            double degen_pressure = 1/8*std::pow(3/consts::PI,1/3)*consts::h*consts::c*std::pow(getDensity()/consts::mH, 4/3);
            pressure = std::max(degen_pressure, pressure);
            if (pressure == degen_pressure) degenerate = 1;
            else degenerate = 0;
        
            return pressure;
        }

        void setPressure(double p) {
            pressure = p;
        }

        void setMass(long double m) {
            if (m > 0) {
                mass = m;
            } else {
                mass = 0;
            }
        }

        void addMass(double m) {
            setMass(mass + m);
        }

        float getX() {
            // hydrogen = hydrogen/(hydrogen+helium+metals);
            return hydrogen;
        }

        void setX(float h) {
            if (h < 0) hydrogen = 0;
            else if (h > 1) hydrogen = 1;
            else hydrogen = h;
            // if (round(mass) > 0) {
                
            // } else {
            //     // std::cout << mass << " " << h << std::endl;
            //     hydrogen = 0;
            //     // hydrogen = h;
            // }

        }

        float getY() {
            // helium = helium/(hydrogen+helium+metals);
            return helium;
        }

        void setY(float he) {
            if (he < 0) helium = 0;
            else if (he > 1) helium = 1;
            else helium = he;
            // if (round(mass) > 0) {
                
            // } else {
            //     helium = 0;
            //     // helium = he;
            // }
        }

        float getZ() {
            // metals = 1 - hydrogen - helium;
            // metals = metals/(hydrogen+helium+metals);
            return metals;
        }

        void setZ(float Z) {
            if (Z < 0) metals = 0;
            else if (Z > 1) metals = 1;
            else metals = Z;
            
        }

        void HtoHe(float dX) {
            if (getX() - dX < 0) {
                setX(0);
                setY(1-getZ());
            } else {
                setX(getX() - dX);
                setY(getY() + dX);
            }
        }

        void HetoZ(float dY) {
            if (getY() - dY < 0) {
                setY(0);
                setZ(1-getX());
            } else {
                setY(getY() - dY);
                setZ(getZ() + dY);
            }
        }

        long double gete(bool recalculate=true) {
            if (recalculate) {
                // if (round(getDensity()*1000) == 0) e = 0;
                // else 
                e = getCv()*getTemp(false)*getDensity();
                // e = E - 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            }
            return e;
        }

        void sete(long double energy) {
            if (mass > 0 && energy > 0) {
                e = energy;
            } else {
                e = 0;
            }
        }
        
        long double getE(bool recalculate=false) {
            // if (recalculate) {
            //     E = e + 1/2*getDensity()*velocity.getMag()*velocity.getMag();
            // }
            return E;
        }

        void setE(long double energy) {
            if (mass > 0 && energy > 0) {
                E = energy;
            } else {
                E = 0;
            }
        }

        uint getDegenerate() {
            return degenerate;
        }

        double getFusionEnergy() {
            return fusionE;
        }

        void setFusionEnergy(double fE) {
            fusionE = fE;
        }

        double getGravPotential() {
            return gravPotential;
        }

        void setGravPotential(double u) {
            gravPotential = u;
        }

        FluidCell *getRight() {
            return neighbors.right;
        }
        FluidCell *getLeft() {
            return neighbors.left;
        }
        FluidCell *getTop() {
            return neighbors.top;
        }
        FluidCell *getBottom() {
            return neighbors.bottom;
        }

        void setRight(FluidCell *neighbor) {
            neighbors.right = neighbor;
        }
        void setLeft(FluidCell *neighbor) {
            neighbors.left = neighbor;
        }
        void setTop(FluidCell *neighbor) {
            neighbors.top = neighbor;
        }
        void setBottom(FluidCell *neighbor) {
            neighbors.bottom = neighbor;
        }

        void transferMass(FluidCell *other, double m) {
            if (other != NULL) { // halt at the boundaries; NOTHING GETS OUT
                if (m > 0) {
                    if (this->getMass() < m) {
                        m = this->getMass();
                    }
                    double otherMass = other->getMass();
                    other->setMass(otherMass + m);
                    this->setMass(mass - m);
                } else if (m < 0) {
                    m = -m;
                    if (other->getMass() < m) {
                        m = other->getMass();
                    }
                    double otherMass = other->getMass();
                    other->setMass(otherMass - m);
                    this->setMass(mass + m);
                }
            }
        }

        void setVelocity(double vx, double vy) {
            velocity.setVx(vx);
            velocity.setVy(vy);
        }

        void setVelocity(VelocityVector v) {
            velocity = v;
        }

        VelocityVector getVelocity() {
            // return VelocityVector(0,1);
            return velocity;
        }

        void refine() {
            refinedThisStep = true;
            if (hasChildren()) return;

            nw = new FluidCell(this, depth+1, (index << 2) | 0, mass/4, width/2, height/2, temperature, hydrogen, helium, metals);
            nw->setGravPotential(gravPotential);
            ne = new FluidCell(this, depth+1, (index << 2) | 1, mass/4, width/2, height/2, temperature, hydrogen, helium, metals);
            ne->setGravPotential(gravPotential);
            sw = new FluidCell(this, depth+1, (index << 2) | 2, mass/4, width/2, height/2, temperature, hydrogen, helium, metals);
            sw->setGravPotential(gravPotential);
            se = new FluidCell(this, depth+1, (index << 2) | 3, mass/4, width/2, height/2, temperature, hydrogen, helium, metals);
            se->setGravPotential(gravPotential);
            nw->sete(e);
            ne->sete(e);
            sw->sete(e);
            se->sete(e);
            nw->setE(E);
            ne->setE(E);
            sw->setE(E);
            se->setE(E);
            nw->setFusionEnergy(fusionE);
            ne->setFusionEnergy(fusionE);
            sw->setFusionEnergy(fusionE);
            se->setFusionEnergy(fusionE);

        }

        // FluidCell **getChildren() {
        //     return children;
        // }

        bool hasChildren() {
            return ((nw != nullptr) || (ne != nullptr) || (sw != nullptr) || (se != nullptr));
        }

        void absorbChildProps() {
            if (isLeaf()) return;
            FluidCell *children[4] {nw, ne, sw, se};
            long double totMass = 0; 
            long double tote = 0;
            long double totE = 0;
            long double totFusE = 0;
            long double avgX = 0;
            long double avgY = 0;
            double avgTemp = 0;
            double avgPressure = 0;
            double avgGravPotential = 0;

            for (int i = 0; i < 4; i++) {
                totMass += children[i]->getMass();
                tote += children[i]->gete(false);
                totE += children[i]->getE(false);
                totFusE += children[i]->getFusionEnergy()*children[i]->getMass();
                avgTemp += children[i]->getTemp(false);
                avgPressure += children[i]->getPressure(false);
                avgGravPotential += children[i]->getGravPotential();
                avgX += children[i]->getX()*children[i]->getMass();
                avgY += children[i]->getY()*children[i]->getMass();
            }

            avgTemp /= 4;
            avgPressure /= 4;
            avgGravPotential /= 4;
            tote /= 4;
            totE /= 4;
            totFusE /= 4;
            avgX /= totMass;
            avgY /= totMass;

            setMass(totMass);
            sete(tote);
            setE(totE);
            setFusionEnergy(totFusE);
            setTemp(avgTemp);
            setPressure(avgPressure);
            setGravPotential(avgGravPotential);
            setX(avgX);
            setY(avgY);
        }

        void coarsen() {
            if (hasChildren()) {
                absorbChildProps();
                // delete nw;
                // delete ne;
                // delete sw;
                // delete se;
                nw = nullptr;
                ne = nullptr;
                sw = nullptr;
                se = nullptr;
            }
        }

        bool isLeaf() {
            return !(hasChildren());
        }

        FluidCell *nw, *ne, *sw, *se;
        FluidCell *parent;
        uint32_t depth, index;
        VelocityVector newVelocity;
        long double newMass, newE;
        long double newMassx, newMassy, newEx, newEy;
        float newX, newY, newZ;
        float newXx, newXy, newYx, newYy, newZx, newZy;
        bool shouldRefine = false;
        bool shouldCoarsen = false;
        bool refinedThisStep = false;
    private:
        long double mass, e, E;
        double size, density, temperature, pressure, gravPotential, width, height, fusionE; 
        uint degenerate;
        int row, col;
        VelocityVector velocity;
        Neighbors neighbors;
        float hydrogen, helium, metals;
        // FluidCell *children[4] = {nw, ne, sw, se};
};

class FluidGrid {
    public:
        FluidGrid(double width, double height, int startDepth, float dt): width(width), height(height), dt(dt) {
            root = new FluidCell(nullptr, 0, 0, 10, width, height, 1, 0.8, 0.15, 0.05);
            int numCells = refineGridtoDepth(root, 0, startDepth);
            assignIDs(root, 0, 0);
            leafCells.reserve(pow(4,maxDepth));
            updateLeafCells();
            setNeighbors();
            
            double mass, temp;
            maxV = 0;
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                // mass = (std::rand() % 100) * 1e16;
                xyPos xy = getXY(cell->depth,cell->index);
                double radius = pow(width*height,0.5)/4;
                double dist = pow(pow(width/2-xy.x,2) + pow(height/2-xy.y,2),0.5);
                double Tc = 4e6;
                // double Tc = 1e8;
                long double mc = 2e21;
                // mass = dist < radius ? pow(1-dist/radius,3)*1e24/numCells : 1e18/numCells;
                mass = dist < radius ? std::sin(consts::PI/2 * dist/radius)/(consts::PI *dist/radius) * mc/numCells : 0/numCells;
                // temp = dist < radius ? pow(1-dist/radius,3)*1e6 : 10;
                // temp = dist < radius ? std::cos(consts::PI/2 * dist / (radius*1.5)) / (consts::PI * dist / (radius*1.5)) * 1e5 : 10;
                temp = dist < radius ? std::sin(consts::PI/2 * dist / radius) / (consts::PI * dist / radius) * Tc : 0;
                // temp = 1e6;
                if (dist == 0) {
                    mass = mc/numCells;
                    temp = Tc;
                }

                cell->setMass(mass);
                cell->setTemp(temp);
                cell->gete(true);
                cell->setE(cell->gete(false));
                cell->getPressure(true);
                if (abs(cell->getVelocity().getVx()) > maxV) maxV = abs(cell->getVelocity().getVx());
                else if (abs(cell->getVelocity().getVy()) > maxV) maxV = abs(cell->getVelocity().getVy());
            }
            
            minSize = std::min(width,height) / pow(pow(4,startDepth),0.5);
            // if (maxV == 0) {
                maxV = 0.4 * minSize/getdt();
            // }
        }

        int refineGridtoDepth(FluidCell *cell, int depth, int mDepth) {
            if (depth >= mDepth) {
                return 1;
            } else {
                cell->refine();
                FluidCell *children[4] {cell->nw, cell->ne, cell->sw, cell->se};
                int numCells = 0;
                for (int i = 0; i < 4; i++) {
                    numCells += refineGridtoDepth(children[i], depth+1, mDepth);
                }
                return numCells;
            }
        }

        FluidCell *adjacentCell(int depth, int index, Direction dir) {
            return adjacentCellHelper(depth, index, dir, depth, index);
        }

        FluidCell *adjacentCell(FluidCell *cell, Direction dir) {
            TreeLoc loc = getLocFromID(getIDfromLeaf(cell));
            int depth = loc.depth;
            int index = loc.index;
            FluidCell *adj = adjacentCellHelper(depth, index, dir, depth, index);
            if (adj == 0) return cell;
            else return adj;
        }

        FluidCell *adjacentCellHelper(int depth, int index, Direction dir, int origDepth, int origIndex) {
            // If depth <= 0, we are at a boundary
            if (depth <= 0) return 0;
            uint64_t id = getID(depth, index);
            uint idx = id & 0x3; // between 0 and 3 (0:nw, 3:se)
            uint nextIdx = -1;
            switch (dir){
                case LEFT:
                    if (idx & 0x1) nextIdx = idx - 1;
                    break;
                case RIGHT:
                    if (!(idx & 0x1)) nextIdx = idx + 1;
                    break;
                case UP:
                    if ((idx >> 1) & 0x1) nextIdx = idx - 2;
                    break;
                case DOWN:
                    if (!((idx >> 1) & 0x1)) nextIdx = idx + 2;
            }

            if (nextIdx != -1) {
                // hide two lowest bits and set them to 0 and then add nextIdx
                uint64_t nextID = getID(depth, ((index >> 2) << 2) | nextIdx);
                Direction opDir;
                FluidCell *currCell;
                FluidCell *cell;
                int vecIdx;// to choose which candidate to use
                // if we are at an internal node
                if (!idToCell.count(nextID)) {
                    
                    switch (dir){
                        case LEFT:
                            opDir = RIGHT;
                            if ((origIndex >> 1) & 0x1) vecIdx = 1;
                            else vecIdx = 0; 
                            break;
                        case RIGHT:
                            opDir = LEFT;
                            if ((origIndex >> 1) & 0x1) vecIdx = 1;
                            else vecIdx = 0;
                            break;
                        case UP:
                            opDir = DOWN;
                            if (origIndex & 0x1) vecIdx = 1;
                            else vecIdx = 0;
                            break;
                        case DOWN:
                            opDir = UP;
                            if (origIndex & 0x1) vecIdx = 1;
                            else vecIdx = 0;
                    }
                    // // Get ID of cell to coarsen
                    // uint64_t parentID = getID(depth, index);
                    // // Index of first (nw) child will be 4 times index of parent
                    // uint32_t childIndex = ((uint32_t) parentID << 2);
                    // // Child is in next layer—get the ID
                    // uint64_t childID = getID(depth+1, childIndex);
                    // // Find child in map and grab it's parent to coarsen
                    // FluidCell *child = idToCell[childID];
                    uint32_t childIdx = ((uint32_t) nextID << 2) & 0xFFFFFFFF;
                    int childDepth = depth+1;
                    uint64_t childID = (((uint64_t) childDepth) << 32) | (uint64_t) childIdx;
                    int nDescended = 1;
                    // while we are not at a leaf, find the ID of a leaf in this
                    //// tree
                    while (!idToCell.count(childID)) {
                        childIdx = ((uint32_t) childIdx << 2) & 0xFFFFFFFF;
                        childDepth += 1;
                        childID = (((uint64_t) childDepth) << 32) | (uint64_t) childIdx;
                        nDescended += 1;
                    }
                    currCell = idToCell[childID];
                    for (int i = 0; i < nDescended; i++) {
                        currCell = currCell->parent;
                    }

                    cell = traverseInDirection(currCell,opDir,depth,origDepth, origIndex);
                    // This would normally happen with getIDfromLeaf but 
                    //// I wanted to hide the unnecessary warning about 
                    //// non-leafage.
                    if (cell->isLeaf()) nextID = getIDfromLeaf(cell); 
                    else nextID = 0;
                }
                uint32_t leafDepth = (nextID >> 32) & 0xFFFFFFFF;
                // if the resolutions are the same or neighbor is more coarse
                if (nextID == 0) {
                    return cell;
                } else if ((leafDepth == origDepth) || (leafDepth == origDepth - 1) || (leafDepth == origDepth + 1)) {
                    return idToCell[nextID];
                    
                // if the neighbor is more refined by a level
                // } else if (leafDepth == origDepth + 1) {
                //     return idToCell[nextID]->parent;
                } else {
                    // std::cout << "dir: " << dir << std::endl;
                    // std::cout << "jumper: ";
                    // printID(nextID);
                    // std::cout << std::endl;
                    // std::cout << origDepth << " ";
                    // printID(origIndex);
                    // std::cout << std::endl;
                    // std::cerr << "Your neighbor has a jump in resolution, and we will not handle that for now. Returning a null pointer.\n";
                    return nullptr;
                }
                
            } else {
                return adjacentCellHelper(depth - 1, index >> 2, dir, origDepth, origIndex);
            }

        }

        FluidCell *traverseInDirection(FluidCell *cell, Direction dir, int depth, int origDepth, int origIdx) {
            /* This will only go as deep as the orignal (target) depth and try to follow the origIdx's location. */
            if (cell->isLeaf() || (depth == origDepth)) {
                return cell;
            }
            uint16_t idx;
            switch (dir){
                case LEFT:
                    if ((origIdx >> (2*(origDepth-depth-1)) >> 1) & 0x1) {
                        return traverseInDirection(cell->sw,dir,depth+1,origDepth,origIdx);
                    } else {
                        return traverseInDirection(cell->nw,dir,depth+1,origDepth,origIdx);
                    }
                    break;
                case RIGHT:
                    if ((origIdx >> (2*(origDepth-depth-1)) >> 1) & 0x1) {
                        return traverseInDirection(cell->se,dir,depth+1,origDepth,origIdx);
                    } else {
                        return traverseInDirection(cell->ne,dir,depth+1,origDepth,origIdx);
                    }
                    break;
                case UP:
                    if (origIdx >> (2*(origDepth-depth-1)) & 0x1) {
                        return traverseInDirection(cell->ne,dir,depth+1,origDepth,origIdx);
                    } else {
                        return traverseInDirection(cell->nw,dir,depth+1,origDepth,origIdx);
                    }
                    break;
                case DOWN:
                    if (origIdx >> (2*(origDepth-depth-1)) & 0x1) {
                        return traverseInDirection(cell->se,dir,depth+1,origDepth,origIdx);
                    } else {
                        return traverseInDirection(cell->sw,dir,depth+1,origDepth,origIdx);
                    }
            }
            return nullptr;
        }

        uint64_t getIDfromLeaf(FluidCell *cell) {
            // for (const auto& pair : idToCell) {
            //     if (pair.second == cell) {
            //         return pair.first;
            //     }
            // }
            return getID(cell->depth, cell->index);
            // std::cerr << "This cell is not a leaf. Returning ID as 0.\n";
            // return 0;
        }

        TreeLoc getLocFromID(uint64_t id) {
            uint32_t index = id & 0xFFFFFFFF;
            uint32_t depth = (id >> 32) & 0xFFFFFFFF;
            TreeLoc loc;
            loc.index = index;
            loc.depth = depth;
            return loc;
        }

        void setNeighbors() {
            uint32_t depth;
            uint32_t index;
            FluidCell *top, *bottom, *left, *right;
            FluidCell *cell;
            for (const auto& pair : idToCell) {
            // for (int i = 0; i < leafCells.size(); i++) {
                // cell = leafCells[i];
                TreeLoc loc = getLocFromID(pair.first);
                depth = loc.depth;
                index = loc.index;
                // depth = cell->depth;
                // index = cell->index;

                top = adjacentCell(depth,index,UP);
                if (top == 0) top = pair.second;
                pair.second->setTop(top);
                
                bottom = adjacentCell(depth,index,DOWN);
                if (bottom == 0) bottom = pair.second;
                pair.second->setBottom(bottom);
                
                left = adjacentCell(depth,index,LEFT);
                if (left == 0) left = pair.second;
                pair.second->setLeft(left);
                
                right = adjacentCell(depth,index,RIGHT);
                if (right == 0) right = pair.second;
                pair.second->setRight(right);
                
            }

        }

        uint64_t getID(int depth, int index) {
            // std::cout << depth << ", " << index << std::endl;
            assert(depth >= 0);
            assert(index < pow(4,depth));
            return ((uint64_t) depth << 32) | index;
        }

        void assignIDs(FluidCell *cell, int depth, int index) {
            uint64_t id = getID(depth, index);
            if (cell != nullptr) {
                FluidCell *children[4] {cell->nw, cell->ne, cell->sw, cell->se};
                if (cell->isLeaf()) idToCell[id] = cell;
                if (cell->hasChildren()) idToCell.erase(id);
                for (int i = 0; i < 4; i++) {
                    assignIDs(children[i], depth+1, (index << 2 | i));
                }
            } else {
                idToCell.erase(id);
            }
        }

        void printID(uint64_t num) {
            if (num > 1) {
                printID(num / 2);
            }
            std::cout << (num % 2);
        }

        void refineCell(int depth, int index) {
            uint64_t id = getID(depth, index);
            assert(idToCell.count(id)); // Can only refine leaf nodes
            FluidCell *cell = idToCell[id];
            cell->refine();
            // make sure any potential jumps get resolved
            if (cell->depth > cell->getRight()->depth && !cell->getRight()->refinedThisStep) refineCell(cell->getRight()->depth, cell->getRight()->index);
            if (cell->depth > cell->getLeft()->depth && !cell->getLeft()->refinedThisStep) refineCell(cell->getLeft()->depth, cell->getLeft()->index);
            if (cell->depth > cell->getTop()->depth && !cell->getTop()->refinedThisStep) refineCell(cell->getTop()->depth, cell->getTop()->index);
            if (cell->depth > cell->getBottom()->depth && !cell->getBottom()->refinedThisStep) refineCell(cell->getBottom()->depth, cell->getBottom()->index);
            // cell->refinedThisStep = true;
            
            // assignIDs(cell, depth, index);
            // this is faster than recursive call
            idToCell.erase(getID(cell->depth, cell->index));
            idToCell[getID(cell->nw->depth, cell->nw->index)] = cell->nw;
            idToCell[getID(cell->ne->depth, cell->ne->index)] = cell->ne;
            idToCell[getID(cell->sw->depth, cell->sw->index)] = cell->sw;
            idToCell[getID(cell->se->depth, cell->se->index)] = cell->se;
            
            if (cell->hasChildren()) {
                if (cell->nw->getWidth() < minSize) minSize = cell->nw->getWidth();
                else if (cell->nw->getHeight() < minSize) minSize = cell->nw->getHeight();
            }
            
            
        }

        void coarsenCell(int depth, int index) {
            // Get ID of cell to coarsen
            uint64_t parentID = getID(depth, index);
            // Index of first (nw) child will be 4 times index of parent
            uint32_t childIndex = ((uint32_t) parentID << 2);
            // Child is in next layer—get the ID
            uint64_t childID = getID(depth+1, childIndex);
            // Find child in map and grab it's parent to coarsen
            FluidCell *child = idToCell[childID];
            child->parent->coarsen();
            assignIDs(child->parent, depth, index);
        }
    
        std::map<uint64_t,FluidCell*> getIDMap() {
            return idToCell;
        }

        
        xyPos getXY(uint32_t depth, uint32_t index) {
            /* Gets position of bottom left corner of cell */
            uint32_t xmask = 0x1;
            uint32_t ymask = 0x2;
            uint16_t xInd = 0x0;
            uint16_t yInd = 0x0;
            for (int i = 0; i < depth; i++) {
                xInd = ((xmask & index) >> i) | xInd;
                xmask = (xmask << 2);
                yInd = ((ymask & index) >> (i+1)) | yInd;
                ymask = (ymask << 2);
            }
            double cellWidth = width / pow(pow(4,depth),0.5);
            double cellHeight = height / pow(pow(4,depth),0.5);
            xyPos xy;
            xy.x = cellWidth*xInd;
            xy.y = height - cellHeight*(yInd+1);
            return xy;
        }

        FluidCell *getRoot() {
            return root;
        }

        float getdt() {
            return dt;
        }

        void solveGravPotential(int iters) {
            int n;
            for (n = 0; n < iters; n++) {
                for (int i = 0; i < leafCells.size(); i++) {
                    double uL, uR, uT, uB = 0;
                    FluidCell *cell = leafCells[i];
                    if ((cell->getTop() == cell) || (cell->getBottom() == cell) || (cell->getLeft() == cell) || (cell->getRight() == cell)) {
                        cell->setGravPotential(0);
                        continue;
                    }
                    double dens = cell->getDensity();
                    double w2 = pow(cell->getWidth(),2);
                    double h2 = pow(cell->getHeight(),2);
                    cell->getTop()->absorbChildProps();
                    cell->getBottom()->absorbChildProps();
                    cell->getLeft()->absorbChildProps();
                    cell->getRight()->absorbChildProps();
                    uT = cell->getTop()->getGravPotential();
                    uB = cell->getBottom()->getGravPotential();
                    uL = cell->getLeft()->getGravPotential();
                    uR = cell->getRight()->getGravPotential();
                    double u = (uL + uR)/(2*(1+w2/h2)) + (uT + uB)/(2*(h2/w2+1)) - 2 * consts::PI * consts::G * dens / (1/w2 + 1/h2);
                    cell->setGravPotential(u);
                }
                
            }
        }

        void updateVelocities() {
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                cell->setE(cell->gete(false) + 0.5*cell->getDensity()*cell->getVelocity().getMag()*cell->getVelocity().getMag());
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellR = cell->getRight();
                VelocityVector newV = VelocityVector(0,0);
                double gravPotT, gravPotB, gravPotL, gravPotR;
                double pressT, pressB, pressL, pressR;
                double distT = (cell->getHeight()/2 + cellT->getHeight()/2);
                double distB = (cell->getHeight()/2 + cellB->getHeight()/2);
                double distL = (cell->getWidth()/2 + cellL->getWidth()/2);
                double distR = (cell->getWidth()/2 + cellR->getWidth()/2);
                if (cellT->hasChildren()) {
                    gravPotT = (cellT->sw->getGravPotential() + cellT->se->getGravPotential())/2;
                    pressT = (cellT->sw->getPressure(true) + cellT->se->getPressure(true))/2;
                    distT = (cell->getHeight()/2 + cellT->sw->getHeight()/2);
                } else {
                    gravPotT = cellT->getGravPotential();
                    pressT = cellT->getPressure(true);
                }

                if (cellB->hasChildren()) {
                    gravPotB = (cellB->nw->getGravPotential() + cellB->ne->getGravPotential())/2;
                    pressB = (cellB->nw->getPressure(true) + cellB->ne->getPressure(true))/2;
                    distB = (cell->getHeight()/2 + cellB->nw->getHeight()/2);
                } else {
                    gravPotB = cellB->getGravPotential();
                    pressB = cellB->getPressure(true);
                }

                if (cellL->hasChildren()) {
                    gravPotL = (cellL->ne->getGravPotential() + cellL->se->getGravPotential())/2;
                    pressL = (cellL->ne->getPressure(true) + cellL->se->getPressure(true))/2;
                    distL = (cell->getWidth()/2 + cellL->ne->getWidth()/2);
                } else {
                    gravPotL = cellL->getGravPotential();
                    pressL = cellL->getPressure(true);
                }

                if (cellR->hasChildren()) {
                    gravPotR = (cellR->nw->getGravPotential() + cellR->sw->getGravPotential())/2;
                    pressR = (cellR->nw->getPressure(true) + cellR->sw->getPressure(true))/2;
                    distR = (cell->getWidth()/2 + cellR->nw->getWidth()/2);
                } else {
                    gravPotR = cellR->getGravPotential();
                    pressR = cellR->getPressure(true);
                }

                double gradUy = (gravPotT - gravPotB)/(distT+distB);
                double gradUx = (gravPotR - gravPotL)/(distL+distR);
                double gradPy = (pressT - pressB)/(distT+distB);
                double gradPx = (pressR - pressL)/(distL+distR);
                VelocityVector v = cell->getVelocity();
                if (cell->getDensity() > 0.00001*maxDensity) {
                    VelocityVector diff = velDiffusion(cell);
                    // newV.setVx(v.getVx() + (-gradUx-gradPx/cell->getDensity())*getdt());
                    newV.setVx(v.getVx() + (-gradUx-gradPx/cell->getDensity() + diff.getVx()/cell->getDensity())*getdt());
                    // newV.setVy(v.getVy() + (-gradUy-gradPy/cell->getDensity())*getdt());
                    newV.setVy(v.getVy() + (-gradUy-gradPy/cell->getDensity() + diff.getVy()/cell->getDensity())*getdt());
                } else {
                    newV.setVx(0);
                    newV.setVy(0);
                }

                cell->setVelocity(newV);
                // if (abs(cell->getVelocity().getVx()) > 1) {
                //     std::cout << cell->depth << ", ";
                //     printID(cell->index);
                //     std::cout << " " << cell->getVelocity().getVx() << ", " << newV.getVx() << " " << -gradPx/cell->getDensity()*getdt() << std::endl;
                // }
                // std::cout << "vel E: " << cell->gete(false) << std::endl;
            }
        }

        VelocityVector velDiffusion(FluidCell *cell) {
            FluidCell *cellT = cell->getTop();
            FluidCell *cellB = cell->getBottom();
            FluidCell *cellL = cell->getLeft();
            FluidCell *cellR = cell->getRight();
            
            double distT = (cell->getHeight()/2 + cellT->getHeight()/2);
            double distB = (cell->getHeight()/2 + cellB->getHeight()/2);
            double distL = (cell->getWidth()/2 + cellL->getWidth()/2);
            double distR = (cell->getWidth()/2 + cellR->getWidth()/2);

            VelocityVector v = cell->getVelocity();
            VelocityVector vR = cellR->getVelocity();
            VelocityVector vL = cellL->getVelocity();
            VelocityVector vB = cellB->getVelocity();
            VelocityVector vT = cellT->getVelocity();

            VelocityVector l,r;

            if (cellT->hasChildren()) {
                // vT = (cellT->sw->getVelocity() + cellT->se->getVelocity())/2;
                distT = (cell->getHeight()/2 + cellT->sw->getHeight()/2);
                l = cellT->sw->getVelocity();
                r = cellT->se->getVelocity();
                vT = (l + r)/2;
            }

            if (cellB->hasChildren()) {
                // vB = (cellB->nw->getVelocity() + cellB->ne->getVelocity())/2;
                distB = (cell->getHeight()/2 + cellB->nw->getHeight()/2);
                l = cellB->nw->getVelocity();
                r = cellB->ne->getVelocity();
                vB = (l + r)/2;
            }

            if (cellL->hasChildren()) {
                // vL = (cellL->ne->getVelocity() + cellL->se->getVelocity())/2;
                distL = (cell->getWidth()/2 + cellL->ne->getWidth()/2);
                l = cellL->ne->getVelocity();
                r = cellL->se->getVelocity();
                vL = (l + r)/2;
            }

            if (cellR->hasChildren()) {
                // vR = (cellR->nw->getVelocity() + cellR->sw->getVelocity())/2;
                distR = (cell->getWidth()/2 + cellR->nw->getWidth()/2);
                l = cellR->nw->getVelocity();
                r = cellR->sw->getVelocity();
                vR = (l + r)/2;
            }

            VelocityVector dvR2 = (vR - v) / distR;
            VelocityVector dvL2 = (v - vL) / distL;
            VelocityVector dvT2 = (vT - v) / distT;
            VelocityVector dvB2 = (v - vB) / distB;
            VelocityVector divX = (dvR2 - dvL2)/(distL+distR);
            VelocityVector divY = (dvT2 - dvB2)/(distT+distB);

            return std::isfinite((divX+divY).getVx()) && std::isfinite((divX+divY).getVy()) ? (divX + divY)*viscosity : VelocityVector(0,0);

        }

        void energyUpdate() {
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                float energyReduction = 1;
                if (cell->getDensity() < 1) {
                    // energyReduction = cell->getDensity()/maxDensity;
                    // energyReduction = 1e-5;
                    // continue;
                }
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellR = cell->getRight();

                double density = cell->getDensity();

                double gravT, gravB, gravL, gravR;
                double pressT, pressB, pressL, pressR;
                VelocityVector vT, vB, vL, vR;
                double distT = (cell->getHeight()/2 + cellT->getHeight()/2);
                double distB = (cell->getHeight()/2 + cellB->getHeight()/2);
                double distL = (cell->getWidth()/2 + cellL->getWidth()/2);
                double distR = (cell->getWidth()/2 + cellR->getWidth()/2);
                VelocityVector l,r;

                if (cellT->hasChildren()) {
                    gravT = (cellT->sw->getGravPotential() + cellT->se->getGravPotential())/2;
                    pressT = (cellT->sw->getPressure(true) + cellT->se->getPressure(true))/2;
                    l = cellT->sw->getVelocity();
                    r = cellT->se->getVelocity();
                    vT = (l + r)/2;
                    distT = (cell->getHeight()/2 + cellT->sw->getHeight()/2);
                } else {
                    gravT = cellT->getGravPotential();
                    pressT = cellT->getPressure(true);
                    vT = cellT->getVelocity();
                }

                if (cellB->hasChildren()) {
                    gravB = (cellB->nw->getGravPotential() + cellB->ne->getGravPotential())/2;
                    pressB = (cellB->nw->getPressure(true) + cellB->ne->getPressure(true))/2;
                    l = cellB->nw->getVelocity();
                    r = cellB->ne->getVelocity();
                    vB = (l + r)/2;
                    distB = (cell->getHeight()/2 + cellB->nw->getHeight()/2);
                } else {
                    gravB = cellB->getGravPotential();
                    pressB = cellB->getPressure(true);
                    vB = cellB->getVelocity();
                }

                if (cellL->hasChildren()) {
                    gravL = (cellL->ne->getGravPotential() + cellL->se->getGravPotential())/2;
                    pressL = (cellL->ne->getPressure(true) + cellL->se->getPressure(true))/2;
                    l = cellL->ne->getVelocity();
                    r = cellL->se->getVelocity();
                    vL = (l + r)/2;
                    distL = (cell->getWidth()/2 + cellL->ne->getWidth()/2);
                } else {
                    gravL = cellL->getGravPotential();
                    pressL = cellL->getPressure(true);
                    vL = cellL->getVelocity();
                }

                if (cellR->hasChildren()) {
                    gravR = (cellR->nw->getGravPotential() + cellR->sw->getGravPotential())/2;
                    pressR = (cellR->nw->getPressure(true) + cellR->sw->getPressure(true))/2;
                    l = cellR->nw->getVelocity();
                    r = cellR->sw->getVelocity();
                    vR = (l + r)/2;
                    distR = (cell->getWidth()/2 + cellR->nw->getWidth()/2);
                } else {
                    gravR = cellR->getGravPotential();
                    pressR = cellR->getPressure(true);
                    vR = cellR->getVelocity();
                }
                long double mass = cell->getMass();

                double pressure = cell->getPressure();
                double E = cell->getE(false);

                // nuclear fusion
                double fusionEnergy = fusionRateH(density, cell->getTemp(false), cell->getX()) * getdt();
                double fusionHeEnergy = fusionRateHe(density, cell->getTemp(false), cell->getY()) * getdt();
                // double fusionHeEnergy = 0;
                cell->setFusionEnergy((fusionEnergy+fusionHeEnergy)/getdt());

                // Update hydrogen and helium abundances
                double dX = std::min((double)1.0,fusionEnergy / (consts::QH_He * density));
                cell->HtoHe(dX);
                double dY = std::min((double)1.0,fusionHeEnergy / (consts::QHe_C * density));
                cell->HetoZ(dY);

                if (dX > 1) fusionEnergy = consts::QH_He * density;
                if (dY > 1) fusionHeEnergy = consts::QHe_C * density;

                E += energyReduction*(fusionEnergy+fusionHeEnergy);

                // pressure-work term
                double div = (vR.getVx()-vL.getVx())/(distR+distL) + (vT.getVy()-vB.getVy())/(distT+distB);
                double pdiv = pressure * div;
                if (E - pdiv * getdt() >= 0) {
                    E += energyReduction*(-pdiv * getdt());
                } else {
                    // std::cerr << "Warning: Pressure work term causing negative energy.\n";
                    E = 0.0f; // Prevent negative energy
                }

                // cooling
                double T = cell->getTemp(false);
                double Z = cell->getZ();
                double lambda = coolingFunction(density,T,Z);
                if (E - lambda * getdt() >= 0) {
                    E += -lambda*getdt();
                } else {
                    E = 0.0f;
                }

                // gravitational-work term
                double grav = cell->getGravPotential();
                // double gravB = getCell(i+1,j)->getGravPotential();
                double gy = -(gravT-gravB)/(distT+distB);
                // double gravR = getCell(i,j+1)->getGravPotential();
                double gx = -(gravR-gravL)/(distR+distL);
                double gravWork = (cell->getVelocity().getVx()*gx+cell->getVelocity().getVy()*gy)*density;
                if (E + gravWork * getdt() >= 0) {
                    E += energyReduction*(gravWork * getdt());
                } else {
                    // std::cerr << "Warning: Gravitational work term causing negative energy.\n";
                    E = 0.0f; // Prevent negative energy
                }
                // std::cout << "E: " << E << std::endl;
                // std::cout << "E: " << cell->getE(false) << std::endl;

                // E += energyReduction*(radiativeTransfer(cell)*getdt());

                cell->setE(E);
                
            }
        }

        long double radiativeTransfer(FluidCell *cell) {
            FluidCell *cellT = cell->getTop();
            FluidCell *cellB = cell->getBottom();
            FluidCell *cellL = cell->getLeft();
            FluidCell *cellR = cell->getRight();

            double density = cell->getDensity();
            double densityT = cellT->getDensity();
            double densityB = cellB->getDensity();
            double densityL = cellL->getDensity();
            double densityR = cellR->getDensity();

            double X = cell->getX();
            double XT = cellT->getX();
            double XB = cellB->getX();
            double XL = cellL->getX();
            double XR = cellR->getX();

            double temp = cell->getTemp(false);
            double tempT = cellT->getTemp(false);
            double tempB = cellB->getTemp(false);
            double tempL = cellL->getTemp(false);
            double tempR = cellR->getTemp(false);

            double E = cell->getE(false);
            double ET = cellT->getE(false);
            double EB = cellB->getE(false);
            double EL = cellL->getE(false);
            double ER = cellR->getE(false);

            double distT = (cell->getHeight()/2 + cellT->getHeight()/2);
            double distB = (cell->getHeight()/2 + cellB->getHeight()/2);
            double distL = (cell->getWidth()/2 + cellL->getWidth()/2);
            double distR = (cell->getWidth()/2 + cellR->getWidth()/2);



            if (cellT->hasChildren()) {
                densityT = (cellT->sw->getDensity() + cellT->se->getDensity())/2;
                tempT = (cellT->sw->getTemp(false) + cellT->se->getTemp(false))/2;
                ET = (cellT->sw->getE(false) + cellT->se->getE(false))/2;
                XT = (cellT->sw->getX() + cellT->se->getX())/2;
                distT = (cell->getHeight()/2 + cellT->sw->getHeight()/2);
            }

            if (cellB->hasChildren()) {
                densityB = (cellB->nw->getDensity() + cellB->ne->getDensity())/2;
                tempB = (cellB->nw->getTemp(false) + cellB->ne->getTemp(false))/2;
                EB = (cellB->nw->getE(false) + cellB->ne->getE(false))/2;
                XB = (cellB->nw->getX() + cellB->ne->getX())/2;
                distB = (cell->getHeight()/2 + cellB->nw->getHeight()/2);
            }

            if (cellL->hasChildren()) {
                densityL = (cellL->ne->getDensity() + cellL->se->getDensity())/2;
                tempL = (cellL->ne->getTemp(false) + cellL->se->getTemp(false))/2;
                EL = (cellL->ne->getE(false) + cellL->se->getE(false))/2;
                XL = (cellL->ne->getX() + cellL->se->getX())/2;
                distL = (cell->getWidth()/2 + cellL->ne->getWidth()/2);
            }

            if (cellR->hasChildren()) {
                densityR = (cellR->nw->getDensity() + cellR->sw->getDensity())/2;
                tempR = (cellR->nw->getTemp(false) + cellR->sw->getTemp(false))/2;
                ER = (cellR->nw->getE(false) + cellR->sw->getE(false))/2;
                XR = (cellR->nw->getX() + cellR->sw->getX())/2;
                distR = (cell->getWidth()/2 + cellR->nw->getWidth()/2);
            } 

            long double D, Dr, Dl, Dt, Db, Dr2, Dl2, Dt2, Db2;
            double kap, kapR, kapL, kapT, kapB, tau, tauR, tauL, tauT, tauB;

            kap = opacity(density,temp,X);
            kapR = opacity(densityR,tempR,XR);
            kapL = opacity(densityL,tempL,XL);
            kapT = opacity(densityT,tempT,XT);
            kapB = opacity(densityB,tempB,XB);
            tau = opticalDepth(kap,density,cell->getWidth());
            tauR = opticalDepth(kapR,densityR,cellR->getWidth());
            tauL = opticalDepth(kapL,densityL,cellL->getWidth());
            tauT = opticalDepth(kapT,densityT,cellT->getHeight());
            tauB = opticalDepth(kapB,densityB,cellB->getHeight());

            // double gradE = std::sqrt( pow((ER-EL)/(distR+distL),2) + pow((ET-EB)/(distT+distB),2) );

            // double R = gradE/(kap*density*E);
            // double lam = 1/R * (1/std::tanh(R) - 1/R);

            // cell-centered D
            D = (4/3*consts::a*consts::c*pow(temp,3)/(kap*density))/(1+pow(tau,-2));
            Dr = (4/3*consts::a*consts::c*pow(tempR,3)/(kapR*densityR))/(1+pow(tauR,-2));
            Dl = (4/3*consts::a*consts::c*pow(tempL,3)/(kapL*densityL))/(1+pow(tauL,-2));
            Dt = (4/3*consts::a*consts::c*pow(tempT,3)/(kapT*densityT))/(1+pow(tauT,-2));
            Db = (4/3*consts::a*consts::c*pow(tempB,3)/(kapB*densityB))/(1+pow(tauB,-2));
            
            // D at bounds with inverse distance-weighted averages
            Dr2 = ((distR - cell->getWidth()/2) * D + cell->getWidth()/2 * Dr) / distR;
            Dl2 = ((distL - cell->getWidth()/2) * D + cell->getWidth()/2 * Dl) / distL;
            Dt2 = ((distT - cell->getHeight()/2) * D + cell->getHeight()/2 * Dt) / distT;
            Db2 = ((distB - cell->getHeight()/2) * D + cell->getHeight()/2 * Db) / distB;

            long double dEx = (Dr2*(tempR-temp)/distR - Dl2*(temp-tempL)/distL)/cell->getWidth();
            long double dEy = (Dt2*(tempT-temp)/distT - Db2*(temp-tempB)/distB)/cell->getHeight();
            long double dE = dEx + dEy;

            if (dE*getdt() > 1e15) {
                printID(getID(cell->depth, cell->index));
                std::cout << " " << opacity(density,temp,X) << " " << density << " " << densityB << " " << D << " " << Dt2 << " " << Db << " " << tauB << "\n T: " << temp << " " << tempT << " " << tempB << " dE:" << dE*getdt() << " dEx,dEy: " << dEx << "," << dEy << std::endl;
            } else {
                // std::cout << opacity(density,temp,X) << " " << density << " " << D << " " << Dt2 << " " << Db << " " << tauB << "\n T: " << temp << " " << tempT << " " << tempB << " dE:" << dE*getdt() << " dEx,dEy: " << dEx << "," << dEy << std::endl;
            }
            
            if (getID(cell->depth, cell->index) == getID(7, 0x0908)) {
                // std::cout << opacity(density,temp,X) << " rho: " << density << " D: " << D << " T: " << temp << " dE: " << dE*getdt() << " dEx,dEy: " << dEx << "," << dEy << std::endl;
            }
            if (!std::isfinite(dE)) {
                // std::cout << opacity(density,temp,X) << " " << density << " " << D << " T: " << temp << " dE:" << dE*getdt() << " " << getdt() << std::endl;
            }

            return dE;

        }

        void advect() {
            bool minSizeSmall = true; 
            long double thisMaxDensity = 0;
            maxV = 0;
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                FluidCell *cellR = cell->getRight();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                VelocityVector vC = cell->getVelocity();
                float X = cell->getX();
                float Y = cell->getY();
                float Z = cell->getZ();
                if (std::max(abs(vC.getVx()),abs(vC.getVy())) > maxV) {
                    maxV = std::max(abs(vC.getVx()),abs(vC.getVy()));
                }
                if (cell->getHeight() == minSize || cell->getWidth() == minSize) minSizeSmall = false;

                long double massFluxR = 0;
                long double massFluxL = 0;
                long double massFluxT = 0;
                long double massFluxB = 0;
                long double EFluxR = 0;
                long double EFluxL = 0;
                long double EFluxT = 0;
                long double EFluxB = 0;
                long double vxFluxR = 0;
                long double vxFluxL = 0;
                long double vyFluxT = 0;
                long double vyFluxB = 0;
                long double XFluxR = 0;
                long double XFluxL = 0;
                long double XFluxT = 0;
                long double XFluxB = 0;
                long double YFluxR = 0;
                long double YFluxL = 0;
                long double YFluxT = 0;
                long double YFluxB = 0;
                long double ZFluxR = 0;
                long double ZFluxL = 0;
                long double ZFluxT = 0;
                long double ZFluxB = 0;
                double vxR;
                double vxL;
                double vyT;
                double vyB;
                long double massR, massL, massT, massB;
                long double ER, EL, ET, EB;
                long double mass = cell->getMass();
                long double E = cell->getE(false);
                float XR, XL, XT, XB, YR, YL, YT, YB, ZR, ZL, ZT, ZB;

                if (cell->getDensity() > thisMaxDensity) thisMaxDensity = cell->getDensity();
                vxR = (cell->getRight()->getVelocity().getVx() + vC.getVx())/2;
                vxL = (cell->getLeft()->getVelocity().getVx() + vC.getVx())/2;
                vyT = (cell->getTop()->getVelocity().getVy() + vC.getVy())/2;
                vyB = (cell->getBottom()->getVelocity().getVy() + vC.getVy())/2;
                XR = cellR->getX();
                XL = cellL->getX();
                XT = cellT->getX();
                XB = cellB->getX();
                YR = cellR->getY();
                YL = cellL->getY();
                YT = cellT->getY();
                YB = cellB->getY();
                ZR = cellR->getZ();
                ZL = cellL->getZ();
                ZT = cellT->getZ();
                ZB = cellB->getZ();

                // not a boundary
                if (cellR != cell) {
                    if (cellR->hasChildren()) {
                        // vxR = 0;
                        // massR = 0;
                        FluidCell *children[2] = {cellR->nw, cellR->sw};
                        for (int i = 0; i < 2; i++) {
                            vxR = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massR = children[i]->getMass();
                            ER = children[i]->getE(false);
                            massFluxR += vxR > 0 ? mass*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            EFluxR += vxR > 0 ? cell->getSize()*E*vxR*getdt()*children[i]->getHeight()/cell->getSize() : children[i]->getSize()*ER*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            vxFluxR += vxR > 0 ? mass*cell->getVelocity().getVx()*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getVelocity().getVx()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            XFluxR += vxR > 0 ? mass*X*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getX()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            YFluxR += vxR > 0 ? mass*Y*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getY()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            ZFluxR += vxR > 0 ? mass*Z*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getZ()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                        }
                        
                    } else {
                        massR = cellR->getMass();
                        ER = cellR->getE(false);
                        massFluxR = vxR > 0 ? mass*vxR*getdt()*cell->getHeight()/cell->getSize() : massR*vxR*getdt()*cell->getHeight()/cellR->getSize();
                        EFluxR = vxR > 0 ? cell->getSize()*E*vxR*getdt()*cell->getHeight()/cell->getSize() : cellR->getSize()*ER*vxR*getdt()*cell->getHeight()/cellR->getSize();
                        vxFluxR = vxR > 0 ? mass*cell->getVelocity().getVx()*vxR*getdt()*cell->getHeight()/cell->getSize() : massR*cellR->getVelocity().getVx()*vxR*getdt()*cell->getHeight()/cellR->getSize();
                        XFluxR = vxR > 0 ? mass*X*vxR*getdt()*cell->getHeight()/cell->getSize() : massR*cellR->getX()*vxR*getdt()*cell->getHeight()/cellR->getSize();
                        YFluxR = vxR > 0 ? mass*Y*vxR*getdt()*cell->getHeight()/cell->getSize() : massR*cellR->getY()*vxR*getdt()*cell->getHeight()/cellR->getSize();
                        ZFluxR = vxR > 0 ? mass*Z*vxR*getdt()*cell->getHeight()/cell->getSize() : massR*cellR->getZ()*vxR*getdt()*cell->getHeight()/cellR->getSize();
                    }
                }
                // not a boundary
                if (cellL != cell) {
                    if (cellL->hasChildren()) {
                        // vxL = 0;
                        // massL = 0;
                        FluidCell *children[2] = {cellL->ne, cellL->se};
                        for (int i = 0; i < 2; i++) {
                            vxL = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massL = children[i]->getMass();
                            EL = children[i]->getE(false);
                            massFluxL += vxL < 0 ? mass*vxL*getdt()*children[i]->getHeight()/cell->getSize() : massL*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                            EFluxL += vxL < 0 ? cell->getSize()*E*vxL*getdt()*children[i]->getHeight()/cell->getSize() : children[i]->getSize()*EL*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                            vxFluxL += vxL < 0 ? mass*cell->getVelocity().getVx()*vxL*getdt()*children[i]->getHeight()/cell->getSize() : massL*children[i]->getVelocity().getVx()*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                            XFluxL += vxL < 0 ? mass*X*vxL*getdt()*children[i]->getHeight()/cell->getSize() : massL*children[i]->getX()*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                            YFluxL += vxL < 0 ? mass*Y*vxL*getdt()*children[i]->getHeight()/cell->getSize() : massL*children[i]->getY()*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                            ZFluxL += vxL < 0 ? mass*Z*vxL*getdt()*children[i]->getHeight()/cell->getSize() : massL*children[i]->getZ()*vxL*getdt()*children[i]->getHeight()/children[i]->getSize();
                        }
                        
                    } else {
                        massL = cellL->getMass();
                        EL = cellL->getE(false);
                        massFluxL = vxL < 0 ? mass*vxL*getdt()*cell->getHeight()/cell->getSize() : massL*vxL*getdt()*cell->getHeight()/cellL->getSize();
                        EFluxL = vxL < 0 ? cell->getSize()*E*vxL*getdt()*cell->getHeight()/cell->getSize() : cellL->getSize()*EL*vxL*getdt()*cell->getHeight()/cellL->getSize();
                        vxFluxL = vxL < 0 ? mass*cell->getVelocity().getVx()*vxL*getdt()*cell->getHeight()/cell->getSize() : massL*cellL->getVelocity().getVx()*vxL*getdt()*cell->getHeight()/cellL->getSize();
                        XFluxL = vxL < 0 ? mass*X*vxL*getdt()*cell->getHeight()/cell->getSize() : massL*cellL->getX()*vxL*getdt()*cell->getHeight()/cellL->getSize();
                        YFluxL = vxL < 0 ? mass*Y*vxL*getdt()*cell->getHeight()/cell->getSize() : massL*cellL->getY()*vxL*getdt()*cell->getHeight()/cellL->getSize();
                        ZFluxL = vxL < 0 ? mass*Z*vxL*getdt()*cell->getHeight()/cell->getSize() : massL*cellL->getZ()*vxL*getdt()*cell->getHeight()/cellL->getSize();
                    }
                }
                // not a boundary
                if (cellT != cell) {
                    if (cellT->hasChildren()) {
                        // vyT = 0;
                        // massT = 0;
                        FluidCell *children[2] = {cellT->sw, cellT->se};
                        for (int i = 0; i < 2; i++) {
                            vyT = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massT = children[i]->getMass();
                            ET = children[i]->getE(false);
                            massFluxT += vyT > 0 ? mass*vyT*getdt()*children[i]->getWidth()/cell->getSize() : massT*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                            EFluxT += vyT > 0 ? cell->getSize()*E*vyT*getdt()*children[i]->getWidth()/cell->getSize() : children[i]->getSize()*ET*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                            vyFluxT += vyT > 0 ? mass*cell->getVelocity().getVy()*vyT*getdt()*children[i]->getWidth()/cell->getSize() : massT*children[i]->getVelocity().getVy()*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                            XFluxT += vyT > 0 ? mass*X*vyT*getdt()*children[i]->getWidth()/cell->getSize() : massT*children[i]->getX()*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                            YFluxT += vyT > 0 ? mass*Y*vyT*getdt()*children[i]->getWidth()/cell->getSize() : massT*children[i]->getY()*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                            ZFluxT += vyT > 0 ? mass*Z*vyT*getdt()*children[i]->getWidth()/cell->getSize() : massT*children[i]->getZ()*vyT*getdt()*children[i]->getWidth()/children[i]->getSize();
                        }
                        
                    } else {
                        massT = cellT->getMass();
                        ET = cellT->getE(false);
                        massFluxT = vyT > 0 ? mass*vyT*getdt()*cell->getWidth()/cell->getSize() : massT*vyT*getdt()*cell->getWidth()/cellT->getSize();
                        EFluxT = vyT > 0 ? cell->getSize()*E*vyT*getdt()*cell->getWidth()/cell->getSize() : cellT->getSize()*ET*vyT*getdt()*cell->getWidth()/cellT->getSize();
                        vyFluxT = vyT > 0 ? mass*cell->getVelocity().getVy()*vyT*getdt()*cell->getWidth()/cell->getSize() : massT*cellT->getVelocity().getVy()*vyT*getdt()*cell->getWidth()/cellT->getSize();
                        XFluxT = vyT > 0 ? mass*X*vyT*getdt()*cell->getWidth()/cell->getSize() : massT*cellT->getX()*vyT*getdt()*cell->getWidth()/cellT->getSize();
                        YFluxT = vyT > 0 ? mass*Y*vyT*getdt()*cell->getWidth()/cell->getSize() : massT*cellT->getY()*vyT*getdt()*cell->getWidth()/cellT->getSize();
                        ZFluxT = vyT > 0 ? mass*Z*vyT*getdt()*cell->getWidth()/cell->getSize() : massT*cellT->getZ()*vyT*getdt()*cell->getWidth()/cellT->getSize();
                    }
                }
                // not a boundary
                if (cellB != cell) {
                    if (cellB->hasChildren()) {
                        // vyB = 0;
                        // massB = 0;
                        FluidCell *children[2] = {cellB->nw, cellB->ne};
                        for (int i = 0; i < 2; i++) {
                            vyB = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massB = children[i]->getMass();
                            EB = children[i]->getE(false);
                            massFluxB += vyB < 0 ? mass*vyB*getdt()*children[i]->getWidth()/cell->getSize() : massB*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                            EFluxB += vyB < 0 ? cell->getSize()*E*vyB*getdt()*children[i]->getWidth()/cell->getSize() : children[i]->getSize()*EB*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                            vyFluxB += vyB < 0 ? mass*cell->getVelocity().getVy()*vyB*getdt()*children[i]->getWidth()/cell->getSize() : massB*children[i]->getVelocity().getVy()*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                            XFluxB += vyB < 0 ? mass*X*vyB*getdt()*children[i]->getWidth()/cell->getSize() : massB*children[i]->getX()*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                            YFluxB += vyB < 0 ? mass*Y*vyB*getdt()*children[i]->getWidth()/cell->getSize() : massB*children[i]->getY()*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                            ZFluxB += vyB < 0 ? mass*Z*vyB*getdt()*children[i]->getWidth()/cell->getSize() : massB*children[i]->getZ()*vyB*getdt()*children[i]->getWidth()/children[i]->getSize();
                        }
                        
                    } else {
                        massB = cellB->getMass();
                        EB = cellB->getE(false);
                        massFluxB = vyB < 0 ? mass*vyB*getdt()*cell->getWidth()/cell->getSize() : massB*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        EFluxB = vyB < 0 ? cell->getSize()*E*vyB*getdt()*cell->getWidth()/cell->getSize() : cellB->getSize()*EB*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        vyFluxB = vyB < 0 ? mass*cell->getVelocity().getVy()*vyB*getdt()*cell->getWidth()/cell->getSize() : massB*cellB->getVelocity().getVy()*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        XFluxB = vyB < 0 ? mass*X*vyB*getdt()*cell->getWidth()/cell->getSize() : massB*cellB->getX()*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        YFluxB = vyB < 0 ? mass*Y*vyB*getdt()*cell->getWidth()/cell->getSize() : massB*cellB->getY()*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        ZFluxB = vyB < 0 ? mass*Z*vyB*getdt()*cell->getWidth()/cell->getSize() : massB*cellB->getZ()*vyB*getdt()*cell->getWidth()/cellB->getSize();
                        // if (!std::isfinite(XFluxB)) {
                        //     std::cout << "nan: " << mass << " " << massB << ", " << massFluxB << ", " << X << ", " << vyB << ", " << cellB->getX() << ", " << getdt() << std::endl;
                        // }
                    }
                }

                if (getID(7,0x2542) == getID(cell->depth, cell->index)) {
                    // std::cout << "mR: " << massFluxR/massR << " mB: " << massFluxB/massB << " vxR: " << vxR << " vyB: " << vyB  << "\n";
                }

                if (!std::isfinite(EFluxR)) {
                    EFluxR = 0;
                    // std::cout << vxR << ", " << massR << std::endl;
                }
                if (!std::isfinite(EFluxL)) {
                    EFluxL = 0;
                    // std::cout << vxL << ", " << massL << std::endl;
                }
                if (!std::isfinite(EFluxT)) {
                    EFluxT = 0;
                    // std::cout << vyT << ", " << massT << std::endl;
                }
                if (!std::isfinite(EFluxB)) {
                    EFluxB = 0;
                    // std::cout << vyB << ", " << massB << std::endl;
                }
                // setAMR(cell);
                cell->newMass = cell->getMass() - massFluxR + massFluxL - massFluxT + massFluxB;
                cell->newE = E + (- EFluxR + EFluxL - EFluxT + EFluxB)/cell->getSize();
                // if (!std::isfinite(cell->newE)) {
                //     std::cout << "nan: " << EFluxR << ", " << EFluxL << ", " << EFluxT << ", " <<  EFluxB << std::endl;
                //     std::cout << cell->depth << " ";
                //     printID(cell->index);
                //     std::cout << std::endl;
                //     cell->newE = cell->getE(false);
                // }
                // std::cout << "E " << E << ", newE " << cell->newE << std::endl;
                double newVx = (mass*cell->getVelocity().getVx() + -vxFluxR + vxFluxL)/cell->newMass;
                if (abs(newVx*getdt()/cell->getWidth()) > 1) {
                    int sign = newVx > 0 ? 1 : -1;
                    newVx = sign*cell->getWidth()/getdt();
                }
                double newVy = (mass*cell->getVelocity().getVy() + -vyFluxT + vyFluxB)/cell->newMass;
                if (abs(newVy*getdt()/cell->getHeight()) > 1) {
                    int sign = newVy > 0 ? 1 : -1;
                    newVy = sign*cell->getHeight()/getdt();
                }
                // if (abs(cell->getVelocity().getVx()) > 1) {
                //     std::cout << cell->depth << ", ";
                //     printID(cell->index);
                //     std::cout << " " << cell->getVelocity().getVx() << ", " << newVx*getdt()/cell->getWidth() << std::endl;
                // }
                cell->newVelocity = VelocityVector(newVx, newVy);
                cell->newX = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxR + XFluxL - XFluxT + XFluxB)/cell->newMass)) );
                cell->newY = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxR + YFluxL - YFluxT + YFluxB)/cell->newMass)) );
                cell->newZ = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxR + ZFluxL - ZFluxT + ZFluxB)/cell->newMass)) );
            }
            maxDensity = thisMaxDensity;
            if (minSizeSmall) minSize *= 2;
            totMass = 0;
            FluidCell *cell;
            for (int i = 0; i < leafCells.size(); i++) {
                cell = leafCells[i];
                cell->setMass(cell->newMass);
                cell->setE(cell->newE);
                totMass += cell->newMass;
                cell->setVelocity(cell->newVelocity);
                cell->setX(std::max(0.0f,std::min(1.0f,cell->newX)));
                cell->setY(std::max(0.0f,std::min(1.0f,cell->newY)));
                cell->setZ(std::max(0.0f,std::min(1.0f,cell->newZ)));
                cell->sete(cell->getE(false) - 0.5*cell->getDensity()*cell->getVelocity().getMag()*cell->getVelocity().getMag());
                cell->getTemp(true);
                cell->getPressure(true);
            }
        }

        void advect2() {
            bool minSizeSmall = true; 
            long double thisMaxDensity = 0;
            maxV = 0;
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                FluidCell *cellR = cell->getRight();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                VelocityVector vC = cell->getVelocity();
                float X = cell->getX();
                float Y = cell->getY();
                float Z = cell->getZ();
                if (std::max(abs(vC.getVx()),abs(vC.getVy())) > maxV) {
                    maxV = std::max(abs(vC.getVx()),abs(vC.getVy()));
                }
                if (cell->getHeight() == minSize || cell->getWidth() == minSize) minSizeSmall = false;

                long double massFluxR = 0;
                long double massFluxL = 0;
                long double massFluxT = 0;
                long double massFluxB = 0;
                long double EFluxR = 0;
                long double EFluxL = 0;
                long double EFluxT = 0;
                long double EFluxB = 0;
                long double vxFluxR = 0;
                long double vxFluxL = 0;
                long double vyFluxT = 0;
                long double vyFluxB = 0;
                long double XFluxR = 0;
                long double XFluxL = 0;
                long double XFluxT = 0;
                long double XFluxB = 0;
                long double YFluxR = 0;
                long double YFluxL = 0;
                long double YFluxT = 0;
                long double YFluxB = 0;
                long double ZFluxR = 0;
                long double ZFluxL = 0;
                long double ZFluxT = 0;
                long double ZFluxB = 0;
                double vxR;
                double vxL;
                double vyT;
                double vyB;
                long double massR, massL, massT, massB;
                long double ER, EL, ET, EB;
                long double mass = cell->getMass();
                long double E = cell->getE(false);
                float XR, XL, XT, XB, YR, YL, YT, YB, ZR, ZL, ZT, ZB;

                if (cell->getDensity() > thisMaxDensity) thisMaxDensity = cell->getDensity();
                vxR = (cell->getRight()->getVelocity().getVx() + vC.getVx())/2;
                vxL = (cell->getLeft()->getVelocity().getVx() + vC.getVx())/2;
                vyT = (cell->getTop()->getVelocity().getVy() + vC.getVy())/2;
                vyB = (cell->getBottom()->getVelocity().getVy() + vC.getVy())/2;
                XR = cellR->getX();
                XL = cellL->getX();
                XT = cellT->getX();
                XB = cellB->getX();
                YR = cellR->getY();
                YL = cellL->getY();
                YT = cellT->getY();
                YB = cellB->getY();
                ZR = cellR->getZ();
                ZL = cellL->getZ();
                ZT = cellT->getZ();
                ZB = cellB->getZ();
                Direction dir;
                FluidCell *nextCell;

                // not a boundary
                if (cellR != cell) {
                    if (cellR->hasChildren()) {
                        // vxR = 0;
                        // massR = 0;
                        FluidCell *children[2] = {cellR->nw, cellR->sw};
                        for (int i = 0; i < 2; i++) {
                            vxR = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massR = children[i]->getMass();
                            ER = children[i]->getE(false);
                            if (vxR < 0) {
                                dir = RIGHT;
                                nextCell = children[i];
                            } else {
                                dir = LEFT;
                                nextCell = cell;
                            }
                            massFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), MASS, dir, 0);
                            EFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), ENERGY, dir, 0);
                            vxFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), HYDROGEN, dir, 0);
                            YFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), HELIUM, dir, 0);
                            ZFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), METALS, dir, 0);
                            // massFluxR += vxR > 0 ? mass*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            // EFluxR += vxR > 0 ? cell->getSize()*E*vxR*getdt()*children[i]->getHeight()/cell->getSize() : children[i]->getSize()*ER*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            // vxFluxR += vxR > 0 ? mass*cell->getVelocity().getVx()*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getVelocity().getVx()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            // XFluxR += vxR > 0 ? mass*X*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getX()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            // YFluxR += vxR > 0 ? mass*Y*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getY()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                            // ZFluxR += vxR > 0 ? mass*Z*vxR*getdt()*children[i]->getHeight()/cell->getSize() : massR*children[i]->getZ()*vxR*getdt()*children[i]->getHeight()/children[i]->getSize();
                        }
                        
                    } else {
                        // massR = cellR->getMass();
                        // ER = cellR->getE(false);
                        if (vxR < 0) {
                            dir = RIGHT;
                            nextCell = cellR;
                        } else {
                            dir = LEFT;
                            nextCell = cell;
                        }
                        massFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), MASS, dir, 0);
                        EFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), ENERGY, dir, 0);
                        vxFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), HYDROGEN, dir, 0);
                        YFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), HELIUM, dir, 0);
                        ZFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), METALS, dir, 0);
                        
                    }
                    // std::cout << massFluxR << std::endl;
                    if (vxR < 0) {
                        massFluxR *= -1;
                        EFluxR *= -1;
                        vxFluxR *= -1;
                        XFluxR *= -1;
                        YFluxR *= -1;
                        ZFluxR *= -1;
                    }
                }

                // not a boundary
                if (cellL != cell) {
                    if (cellL->hasChildren()) {
                        // vxL = 0;
                        // massL = 0;
                        FluidCell *children[2] = {cellL->ne, cellL->se};
                        for (int i = 0; i < 2; i++) {
                            vxL = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massL = children[i]->getMass();
                            EL = children[i]->getE(false);
                            if (vxL < 0) {
                                dir = RIGHT;
                                nextCell = cell;
                            } else {
                                dir = LEFT;
                                nextCell = children[i];
                            }
                            massFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), MASS, dir, 0);
                            EFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), ENERGY, dir, 0);
                            vxFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), HYDROGEN, dir, 0);
                            YFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), HELIUM, dir, 0);
                            ZFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), METALS, dir, 0);
                        }
                    } else {
                        massL = cellL->getMass();
                        EL = cellL->getE(false);
                        if (vxL < 0) {
                            dir = RIGHT;
                            nextCell = cell;
                        } else {
                            dir = LEFT;
                            nextCell = cellL;
                        }
                        massFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), MASS, dir, 0);
                        EFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), ENERGY, dir, 0);
                        vxFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), HYDROGEN, dir, 0);
                        YFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), HELIUM, dir, 0);
                        ZFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), METALS, dir, 0);
                    }
                    if (vxL < 0) {
                        massFluxL *= -1;
                        EFluxL *= -1;
                        vxFluxL *= -1;
                        XFluxL *= -1;
                        YFluxL *= -1;
                        ZFluxL *= -1;
                    }
                }

                // not a boundary
                if (cellT != cell) {
                    if (cellT->hasChildren()) {
                        // vyT = 0;
                        // massT = 0;
                        FluidCell *children[2] = {cellT->sw, cellT->se};
                        for (int i = 0; i < 2; i++) {
                            vyT = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massT = children[i]->getMass();
                            ET = children[i]->getE(false);
                            if (vyT < 0) {
                                dir = UP;
                                nextCell = children[i];
                            } else {
                                dir = DOWN;
                                nextCell = cell;
                            }
                            massFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), MASS, dir, 0);
                            EFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), ENERGY, dir, 0);
                            vyFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), HYDROGEN, dir, 0);
                            YFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), HELIUM, dir, 0);
                            ZFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), METALS, dir, 0);
                        }
                        
                    } else {
                        massT = cellT->getMass();
                        ET = cellT->getE(false);
                        if (vyT < 0) {
                            dir = UP;
                            nextCell = cellT;
                        } else {
                            dir = DOWN;
                            nextCell = cell;
                        }
                        massFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), MASS, dir, 0);
                        EFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), ENERGY, dir, 0);
                        vyFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), HYDROGEN, dir, 0);
                        YFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), HELIUM, dir, 0);
                        ZFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), METALS, dir, 0);
                    }
                    if (vyT < 0) {
                        massFluxT *= -1;
                        EFluxT *= -1;
                        vyFluxT *= -1;
                        XFluxT *= -1;
                        YFluxT *= -1;
                        ZFluxT *= -1;
                    }
                }
                // not a boundary
                if (cellB != cell) {
                    if (cellB->hasChildren()) {
                        // vyB = 0;
                        // massB = 0;
                        FluidCell *children[2] = {cellB->nw, cellB->ne};
                        for (int i = 0; i < 2; i++) {
                            vyB = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massB = children[i]->getMass();
                            EB = children[i]->getE(false);
                            if (vyB < 0) {
                                dir = UP;
                                nextCell = cell;
                            } else {
                                dir = DOWN;
                                nextCell = children[i];
                            }
                            massFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), MASS, dir, 0);
                            EFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), ENERGY, dir, 0);
                            vyFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), HYDROGEN, dir, 0);
                            YFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), HELIUM, dir, 0);
                            ZFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), METALS, dir, 0);
                        }
                        
                    } else {
                        massB = cellB->getMass();
                        EB = cellB->getE(false);
                        if (vyB < 0) {
                            dir = UP;
                            nextCell = cell;
                        } else {
                            dir = DOWN;
                            nextCell = cellB;
                        }
                        massFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), MASS, dir, 0);
                        EFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), ENERGY, dir, 0);
                        vyFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), HYDROGEN, dir, 0);
                        YFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), HELIUM, dir, 0);
                        ZFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), METALS, dir, 0);
                    }
                    if (vyB < 0) {
                        massFluxB *= -1;
                        EFluxB *= -1;
                        vyFluxB *= -1;
                        XFluxB *= -1;
                        YFluxB *= -1;
                        ZFluxB *= -1;
                    }
                }

                if (!std::isfinite(EFluxR)) {
                    EFluxR = 0;
                    // std::cout << vxR << ", " << massR << std::endl;
                }
                if (!std::isfinite(EFluxL)) {
                    EFluxL = 0;
                    // std::cout << vxL << ", " << massL << std::endl;
                }
                if (!std::isfinite(EFluxT)) {
                    EFluxT = 0;
                    // std::cout << vyT << ", " << massT << std::endl;
                }
                if (!std::isfinite(EFluxB)) {
                    EFluxB = 0;
                    // std::cout << vyB << ", " << massB << std::endl;
                }
                if (getID(cell->depth, cell->index) == getID(7, 0x3094)) {
                    // std::cout << massFluxR << " " << massFluxL << " " << massFluxT << " " << massFluxB << " " << vxR*getdt()/cell->getWidth() << "\n";

                }
                cell->newMass = cell->getMass() - massFluxR + massFluxL - massFluxT + massFluxB;
                cell->newE = E + (- EFluxR + EFluxL - EFluxT + EFluxB)/cell->getSize();
            
                double newVx = (mass*cell->getVelocity().getVx() + -vxFluxR + vxFluxL)/cell->newMass;
                double newVy = (mass*cell->getVelocity().getVy() + -vyFluxT + vyFluxB)/cell->newMass;
                
                cell->newVelocity = VelocityVector(newVx, newVy);
                cell->newX = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxR + XFluxL - XFluxT + XFluxB)/cell->newMass)) );
                cell->newY = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxR + YFluxL - YFluxT + YFluxB)/cell->newMass)) );
                cell->newZ = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxR + ZFluxL - ZFluxT + ZFluxB)/cell->newMass)) );
            }
            maxDensity = thisMaxDensity;
            if (minSizeSmall) minSize *= 2;
            totMass = 0;
            FluidCell *cell;
            for (int i = 0; i < leafCells.size(); i++) {
                cell = leafCells[i];
                if (cell->newMass/cell->getSize() > 100*cell->getDensity()) {
                    // std::cout << cell->depth << " ";
                    // printID(cell->index);
                    // std::cout << std::endl;
                }
                cell->setMass(cell->newMass);
                cell->setE(cell->newE);
                totMass += cell->newMass;
                cell->setVelocity(cell->newVelocity);
                cell->setX(std::max(0.0f,std::min(1.0f,cell->newX)));
                cell->setY(std::max(0.0f,std::min(1.0f,cell->newY)));
                cell->setZ(std::max(0.0f,std::min(1.0f,cell->newZ)));
                cell->sete(cell->getE(false) - 0.5*cell->getDensity()*cell->getVelocity().getMag()*cell->getVelocity().getMag());
                cell->getTemp(true);
                cell->getPressure(true);
            }
        }

        long double getFlux(FluidCell *cell, long double dist, long double boundary, Property prop, Direction dir, long double totFlux) {
            if (!std::isfinite(dist)) return 0;
            long double leng = (dir == LEFT) || (dir == RIGHT) ? cell->getWidth() : cell->getHeight();
            if (cell->hasChildren()) {
                switch (dir) {
                    case RIGHT:
                        return getFlux(cell->nw, dist, boundary, prop, dir, totFlux) + getFlux(cell->sw, dist, boundary, prop, dir, totFlux);
                        break;
                    case LEFT:
                        return getFlux(cell->ne, dist, boundary, prop, dir, totFlux) + getFlux(cell->se, dist, boundary, prop, dir, totFlux);
                        break;
                    case UP:
                        return getFlux(cell->se, dist, boundary, prop, dir, totFlux) + getFlux(cell->sw, dist, boundary, prop, dir, totFlux);
                        break;
                    case DOWN:
                        return getFlux(cell->ne, dist, boundary, prop, dir, totFlux) + getFlux(cell->nw, dist, boundary, prop, dir, totFlux);
                }
            } else if (dist < leng) {
                switch (prop) {
                    case MASS:
                        return totFlux + cell->getMass()*boundary*dist/cell->getSize();
                        break;
                    case MOMENTUM:
                        if (dir == LEFT || dir == RIGHT) {
                            return totFlux + cell->getVelocity().getVx()*cell->getMass()*boundary*dist/cell->getSize();
                        } else {
                            return totFlux + cell->getVelocity().getVy()*cell->getMass()*boundary*dist/cell->getSize();
                        }
                        break;
                    case ENERGY:
                        return cell->getSize()*cell->getE(false)*boundary*dist/cell->getSize();
                        break;
                    case HYDROGEN:
                        return totFlux + cell->getX()*cell->getMass()*boundary*dist/cell->getSize();
                        break;
                    case HELIUM:
                        return totFlux + cell->getY()*cell->getMass()*boundary*dist/cell->getSize();
                        break;
                    case METALS:
                        return totFlux + cell->getZ()*cell->getMass()*boundary*dist/cell->getSize();
                }
            } else {
                FluidCell *nextCell;
                switch (dir) {
                    case RIGHT:
                        nextCell = cell->getRight();
                        break;
                    case LEFT:
                        nextCell = cell->getLeft();
                        break;
                    case UP:
                        nextCell = cell->getTop();
                        break;
                    case DOWN:
                        nextCell = cell->getBottom();
                }
                switch (prop) {
                    case MASS:
                        return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass());
                        break;
                    case MOMENTUM:
                        if (dir == LEFT || dir == RIGHT) {
                            return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass()*cell->getVelocity().getVx());
                        } else {
                            return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass()*cell->getVelocity().getVy());
                        }
                        break;
                    case ENERGY:
                        return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getE(false)*cell->getSize());
                        break;
                    case HYDROGEN:
                        return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass()*cell->getX());
                        break;
                    case HELIUM:
                        return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass()*cell->getY());
                        break;
                    case METALS:
                        return getFlux(nextCell, dist-leng, boundary, prop, dir, totFlux + cell->getMass()*cell->getZ());
                }
            }
            return 0;
        }

        void advect3() {
            bool minSizeSmall = true; 
            long double thisMaxDensity = 0;
            maxV = 0;
            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                FluidCell *cellR = cell->getRight();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                VelocityVector vC = cell->getVelocity();
                float X = cell->getX();
                float Y = cell->getY();
                float Z = cell->getZ();
                if (std::max(abs(vC.getVx()),abs(vC.getVy())) > maxV) {
                    maxV = std::max(abs(vC.getVx()),abs(vC.getVy()));
                }
                if (cell->getHeight() == minSize || cell->getWidth() == minSize) minSizeSmall = false;

                long double massFluxR = 0;
                long double massFluxL = 0;
                long double massFluxT = 0;
                long double massFluxB = 0;
                long double EFluxR = 0;
                long double EFluxL = 0;
                long double EFluxT = 0;
                long double EFluxB = 0;
                long double vxFluxR = 0;
                long double vxFluxL = 0;
                long double vyFluxT = 0;
                long double vyFluxB = 0;
                long double XFluxR = 0;
                long double XFluxL = 0;
                long double XFluxT = 0;
                long double XFluxB = 0;
                long double YFluxR = 0;
                long double YFluxL = 0;
                long double YFluxT = 0;
                long double YFluxB = 0;
                long double ZFluxR = 0;
                long double ZFluxL = 0;
                long double ZFluxT = 0;
                long double ZFluxB = 0;
                double vxR;
                double vxL;
                double vyT;
                double vyB;
                long double massR, massL, massT, massB;
                long double ER, EL, ET, EB;
                long double mass = cell->getMass();
                long double E = cell->getE(false);
                float XR, XL, XT, XB, YR, YL, YT, YB, ZR, ZL, ZT, ZB;

                if (cell->getDensity() > thisMaxDensity) thisMaxDensity = cell->getDensity();
                vxR = (cell->getRight()->getVelocity().getVx() + vC.getVx())/2;
                vxL = (cell->getLeft()->getVelocity().getVx() + vC.getVx())/2;
                vyT = (cell->getTop()->getVelocity().getVy() + vC.getVy())/2;
                vyB = (cell->getBottom()->getVelocity().getVy() + vC.getVy())/2;
                XR = cellR->getX();
                XL = cellL->getX();
                XT = cellT->getX();
                XB = cellB->getX();
                YR = cellR->getY();
                YL = cellL->getY();
                YT = cellT->getY();
                YB = cellB->getY();
                ZR = cellR->getZ();
                ZL = cellL->getZ();
                ZT = cellT->getZ();
                ZB = cellB->getZ();
                Direction dir;
                FluidCell *nextCell;

                // not a boundary
                if (cellR != cell) {
                    if (cellR->hasChildren()) {
                        // vxR = 0;
                        // massR = 0;
                        FluidCell *children[2] = {cellR->nw, cellR->sw};
                        for (int i = 0; i < 2; i++) {
                            vxR = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massR = children[i]->getMass();
                            ER = children[i]->getE(false);
                            if (vxR < 0) {
                                dir = RIGHT;
                                nextCell = children[i];
                            } else {
                                dir = LEFT;
                                nextCell = cell;
                            }
                            massFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                            EFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                            vxFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                            YFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                            ZFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        // massR = cellR->getMass();
                        // ER = cellR->getE(false);
                        if (vxR < 0) {
                            dir = RIGHT;
                            nextCell = cellR;
                        } else {
                            dir = LEFT;
                            nextCell = cell;
                        }
                        massFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                        EFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                        vxFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                        YFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                        ZFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                        
                    }
                    // std::cout << massFluxR << std::endl;
                    if (vxR < 0) {
                        massFluxR *= -1;
                        EFluxR *= -1;
                        vxFluxR *= -1;
                        XFluxR *= -1;
                        YFluxR *= -1;
                        ZFluxR *= -1;
                    }
                }

                // not a boundary
                if (cellL != cell) {
                    if (cellL->hasChildren()) {
                        // vxL = 0;
                        // massL = 0;
                        FluidCell *children[2] = {cellL->ne, cellL->se};
                        for (int i = 0; i < 2; i++) {
                            vxL = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massL = children[i]->getMass();
                            EL = children[i]->getE(false);
                            if (vxL < 0) {
                                dir = RIGHT;
                                nextCell = cell;
                            } else {
                                dir = LEFT;
                                nextCell = children[i];
                            }
                            massFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                            EFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                            vxFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                            YFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                            ZFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                            
                        }
                    } else {
                        massL = cellL->getMass();
                        EL = cellL->getE(false);
                        if (vxL < 0) {
                            dir = RIGHT;
                            nextCell = cell;
                        } else {
                            dir = LEFT;
                            nextCell = cellL;
                        }
                        massFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                        EFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                        vxFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                        YFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                        ZFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                    }
                    if (vxL < 0) {
                        massFluxL *= -1;
                        EFluxL *= -1;
                        vxFluxL *= -1;
                        XFluxL *= -1;
                        YFluxL *= -1;
                        ZFluxL *= -1;
                    }
                }

                // not a boundary
                if (cellT != cell) {
                    if (cellT->hasChildren()) {
                        // vyT = 0;
                        // massT = 0;
                        FluidCell *children[2] = {cellT->sw, cellT->se};
                        for (int i = 0; i < 2; i++) {
                            vyT = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massT = children[i]->getMass();
                            ET = children[i]->getE(false);
                            if (vyT < 0) {
                                dir = UP;
                                nextCell = children[i];
                            } else {
                                dir = DOWN;
                                nextCell = cell;
                            }
                            massFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                            EFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                            vyFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                            YFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                            ZFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        massT = cellT->getMass();
                        ET = cellT->getE(false);
                        if (vyT < 0) {
                            dir = UP;
                            nextCell = cellT;
                        } else {
                            dir = DOWN;
                            nextCell = cell;
                        }
                        massFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                        EFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                        vyFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                        YFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                        ZFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                    }
                    if (vyT < 0) {
                        massFluxT *= -1;
                        EFluxT *= -1;
                        vyFluxT *= -1;
                        XFluxT *= -1;
                        YFluxT *= -1;
                        ZFluxT *= -1;
                    }
                }
                // not a boundary
                if (cellB != cell) {
                    if (cellB->hasChildren()) {
                        // vyB = 0;
                        // massB = 0;
                        FluidCell *children[2] = {cellB->nw, cellB->ne};
                        for (int i = 0; i < 2; i++) {
                            vyB = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massB = children[i]->getMass();
                            EB = children[i]->getE(false);
                            if (vyB < 0) {
                                dir = UP;
                                nextCell = cell;
                            } else {
                                dir = DOWN;
                                nextCell = children[i];
                            }
                            massFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                            EFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                            vyFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                            YFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                            ZFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        massB = cellB->getMass();
                        EB = cellB->getE(false);
                        if (vyB < 0) {
                            dir = UP;
                            nextCell = cell;
                        } else {
                            dir = DOWN;
                            nextCell = cellB;
                        }
                        massFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getMass();}, dir, 0);
                        EFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getE(false)*c->getSize();}, dir, 0);
                        vyFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getX()*c->getMass();}, dir, 0);
                        YFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getY()*c->getMass();}, dir, 0);
                        ZFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c) {return c->getZ()*c->getMass();}, dir, 0);
                    }
                    if (vyB < 0) {
                        massFluxB *= -1;
                        EFluxB *= -1;
                        vyFluxB *= -1;
                        XFluxB *= -1;
                        YFluxB *= -1;
                        ZFluxB *= -1;
                    }
                }

                if (!std::isfinite(EFluxR)) {
                    EFluxR = 0;
                    // std::cout << vxR << ", " << massR << std::endl;
                }
                if (!std::isfinite(EFluxL)) {
                    EFluxL = 0;
                    // std::cout << vxL << ", " << massL << std::endl;
                }
                if (!std::isfinite(EFluxT)) {
                    EFluxT = 0;
                    // std::cout << vyT << ", " << massT << std::endl;
                }
                if (!std::isfinite(EFluxB)) {
                    EFluxB = 0;
                    // std::cout << vyB << ", " << massB << std::endl;
                }

                cell->newMass = cell->getMass() - massFluxR + massFluxL - massFluxT + massFluxB;
                cell->newMassx = cell->getMass() - massFluxR + massFluxL;
                cell->newMassy = cell->getMass() - massFluxT + massFluxB;

                cell->newE = E + (-EFluxR + EFluxL - EFluxT + EFluxB)/cell->getSize();
                cell->newEx = E + (-EFluxR + EFluxL)/cell->getSize();
                cell->newEy = E + (-EFluxT + EFluxB)/cell->getSize();

                double newVx = (mass*cell->getVelocity().getVx() + -vxFluxR + vxFluxL)/cell->newMass;
                double newVy = (mass*cell->getVelocity().getVy() + -vyFluxT + vyFluxB)/cell->newMass;
                
                cell->newVelocity = VelocityVector(newVx, newVy);
                cell->newX = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxR + XFluxL - XFluxT + XFluxB)/cell->newMass)) );
                cell->newXx = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxR + XFluxL)/cell->newMass)) );
                cell->newXy = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxT + XFluxB)/cell->newMass)) );

                cell->newY = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxR + YFluxL - YFluxT + YFluxB)/cell->newMass)) );
                cell->newYx = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxR + YFluxL)/cell->newMass)) );
                cell->newYy = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxT + YFluxB)/cell->newMass)) );

                cell->newZ = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxR + ZFluxL - ZFluxT + ZFluxB)/cell->newMass)) );
                cell->newZx = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxR + ZFluxL)/cell->newMass)) );
                cell->newZy = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxT + ZFluxB)/cell->newMass)) );

            }
            maxDensity = thisMaxDensity;
            if (minSizeSmall) minSize *= 2;

            for (int i = 0; i < leafCells.size(); i++) {
                FluidCell *cell = leafCells[i];
                FluidCell *cellR = cell->getRight();
                FluidCell *cellL = cell->getLeft();
                FluidCell *cellT = cell->getTop();
                FluidCell *cellB = cell->getBottom();
                VelocityVector vC = cell->getVelocity();
                float X = cell->getX();
                float Y = cell->getY();
                float Z = cell->getZ();
                if (std::max(abs(vC.getVx()),abs(vC.getVy())) > maxV) {
                    maxV = std::max(abs(vC.getVx()),abs(vC.getVy()));
                }
                if (cell->getHeight() == minSize || cell->getWidth() == minSize) minSizeSmall = false;

                long double massFluxR = 0;
                long double massFluxL = 0;
                long double massFluxT = 0;
                long double massFluxB = 0;
                long double EFluxR = 0;
                long double EFluxL = 0;
                long double EFluxT = 0;
                long double EFluxB = 0;
                long double vxFluxR = 0;
                long double vxFluxL = 0;
                long double vyFluxT = 0;
                long double vyFluxB = 0;
                long double XFluxR = 0;
                long double XFluxL = 0;
                long double XFluxT = 0;
                long double XFluxB = 0;
                long double YFluxR = 0;
                long double YFluxL = 0;
                long double YFluxT = 0;
                long double YFluxB = 0;
                long double ZFluxR = 0;
                long double ZFluxL = 0;
                long double ZFluxT = 0;
                long double ZFluxB = 0;
                double vxR;
                double vxL;
                double vyT;
                double vyB;
                long double massR, massL, massT, massB;
                long double ER, EL, ET, EB;
                long double mass = cell->getMass();
                long double E = cell->getE(false);
                float XR, XL, XT, XB, YR, YL, YT, YB, ZR, ZL, ZT, ZB;

                if (cell->getDensity() > thisMaxDensity) thisMaxDensity = cell->getDensity();
                vxR = (cell->getRight()->getVelocity().getVx() + vC.getVx())/2;
                vxL = (cell->getLeft()->getVelocity().getVx() + vC.getVx())/2;
                vyT = (cell->getTop()->getVelocity().getVy() + vC.getVy())/2;
                vyB = (cell->getBottom()->getVelocity().getVy() + vC.getVy())/2;
                XR = cellR->getX();
                XL = cellL->getX();
                XT = cellT->getX();
                XB = cellB->getX();
                YR = cellR->getY();
                YL = cellL->getY();
                YT = cellT->getY();
                YB = cellB->getY();
                ZR = cellR->getZ();
                ZL = cellL->getZ();
                ZT = cellT->getZ();
                ZB = cellB->getZ();
                Direction dir;
                FluidCell *nextCell;

                // not a boundary
                if (cellR != cell) {
                    if (cellR->hasChildren()) {
                        // vxR = 0;
                        // massR = 0;
                        FluidCell *children[2] = {cellR->nw, cellR->sw};
                        for (int i = 0; i < 2; i++) {
                            vxR = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massR = children[i]->getMass();
                            ER = children[i]->getE(false);
                            if (vxR < 0) {
                                dir = RIGHT;
                                nextCell = children[i];
                            } else {
                                dir = LEFT;
                                nextCell = cell;
                            }
                            massFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newMassy;}, dir, 0);
                            EFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newEy*c->getSize();}, dir, 0);
                            vxFluxR += getFlux(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newXy*c->getMass();}, dir, 0);
                            YFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newYy*c->getMass();}, dir, 0);
                            ZFluxR += getFlux2(nextCell, std::abs(vxR*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newZy*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        // massR = cellR->getMass();
                        // ER = cellR->getE(false);
                        if (vxR < 0) {
                            dir = RIGHT;
                            nextCell = cellR;
                        } else {
                            dir = LEFT;
                            nextCell = cell;
                        }
                        massFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c){return c->newMassy;}, dir, 0);
                        EFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c){return c->newEy*c->getSize();}, dir, 0);
                        vxFluxR = getFlux(nextCell, std::abs(vxR*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c){return c->newXy*c->getMass();}, dir, 0);
                        YFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c){return c->newYy*c->getMass();}, dir, 0);
                        ZFluxR = getFlux2(nextCell, std::abs(vxR*getdt()), cell->getHeight(), [](FluidCell *c){return c->newZy*c->getMass();}, dir, 0);
                        
                    }
                    // std::cout << massFluxR << std::endl;
                    if (vxR < 0) {
                        massFluxR *= -1;
                        EFluxR *= -1;
                        vxFluxR *= -1;
                        XFluxR *= -1;
                        YFluxR *= -1;
                        ZFluxR *= -1;
                    }
                }

                // not a boundary
                if (cellL != cell) {
                    if (cellL->hasChildren()) {
                        // vxL = 0;
                        // massL = 0;
                        FluidCell *children[2] = {cellL->ne, cellL->se};
                        for (int i = 0; i < 2; i++) {
                            vxL = (children[i]->getVelocity().getVx() + vC.getVx())/2;
                            massL = children[i]->getMass();
                            EL = children[i]->getE(false);
                            if (vxL < 0) {
                                dir = RIGHT;
                                nextCell = cell;
                            } else {
                                dir = LEFT;
                                nextCell = children[i];
                            }
                            massFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newMassy;}, dir, 0);
                            EFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newEy*c->getSize();}, dir, 0);
                            vxFluxL += getFlux(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), MOMENTUM, dir, 0);
                            XFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newXy*c->getMass();}, dir, 0);
                            YFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newYy*c->getMass();}, dir, 0);
                            ZFluxL += getFlux2(nextCell, std::abs(vxL*getdt()), children[i]->getHeight(), [](FluidCell *c){return c->newZy*c->getMass();}, dir, 0);
                            
                        }
                    } else {
                        massL = cellL->getMass();
                        EL = cellL->getE(false);
                        if (vxL < 0) {
                            dir = RIGHT;
                            nextCell = cell;
                        } else {
                            dir = LEFT;
                            nextCell = cellL;
                        }
                        massFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c){return c->newMassy;}, dir, 0);
                        EFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c){return c->newEy*c->getSize();}, dir, 0);
                        vxFluxL = getFlux(nextCell, std::abs(vxL*getdt()), cell->getHeight(), MOMENTUM, dir, 0);
                        XFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c){return c->newXy*c->getMass();}, dir, 0);
                        YFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c){return c->newYy*c->getMass();}, dir, 0);
                        ZFluxL = getFlux2(nextCell, std::abs(vxL*getdt()), cell->getHeight(), [](FluidCell *c){return c->newZy*c->getMass();}, dir, 0);
                    }
                    if (vxL < 0) {
                        massFluxL *= -1;
                        EFluxL *= -1;
                        vxFluxL *= -1;
                        XFluxL *= -1;
                        YFluxL *= -1;
                        ZFluxL *= -1;
                    }
                }

                // not a boundary
                if (cellT != cell) {
                    if (cellT->hasChildren()) {
                        // vyT = 0;
                        // massT = 0;
                        FluidCell *children[2] = {cellT->sw, cellT->se};
                        for (int i = 0; i < 2; i++) {
                            vyT = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massT = children[i]->getMass();
                            ET = children[i]->getE(false);
                            if (vyT < 0) {
                                dir = UP;
                                nextCell = children[i];
                            } else {
                                dir = DOWN;
                                nextCell = cell;
                            }
                            massFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newMassx;}, dir, 0);
                            EFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newEx*c->getSize();}, dir, 0);
                            vyFluxT += getFlux(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newXx*c->getMass();}, dir, 0);
                            YFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newYx*c->getMass();}, dir, 0);
                            ZFluxT += getFlux2(nextCell, std::abs(vyT*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newZx*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        massT = cellT->getMass();
                        ET = cellT->getE(false);
                        if (vyT < 0) {
                            dir = UP;
                            nextCell = cellT;
                        } else {
                            dir = DOWN;
                            nextCell = cell;
                        }
                        massFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c){return c->newMassx;}, dir, 0);
                        EFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c){return c->newEx*c->getSize();}, dir, 0);
                        vyFluxT = getFlux(nextCell, std::abs(vyT*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c){return c->newXx*c->getMass();}, dir, 0);
                        YFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c){return c->newYx*c->getMass();}, dir, 0);
                        ZFluxT = getFlux2(nextCell, std::abs(vyT*getdt()), cell->getWidth(), [](FluidCell *c){return c->newZx*c->getMass();}, dir, 0);
                    }
                    if (vyT < 0) {
                        massFluxT *= -1;
                        EFluxT *= -1;
                        vyFluxT *= -1;
                        XFluxT *= -1;
                        YFluxT *= -1;
                        ZFluxT *= -1;
                    }
                }
                // not a boundary
                if (cellB != cell) {
                    if (cellB->hasChildren()) {
                        // vyB = 0;
                        // massB = 0;
                        FluidCell *children[2] = {cellB->nw, cellB->ne};
                        for (int i = 0; i < 2; i++) {
                            vyB = (children[i]->getVelocity().getVy() + vC.getVy())/2;
                            massB = children[i]->getMass();
                            EB = children[i]->getE(false);
                            if (vyB < 0) {
                                dir = UP;
                                nextCell = cell;
                            } else {
                                dir = DOWN;
                                nextCell = children[i];
                            }
                            massFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newMassx;}, dir, 0);
                            EFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newEx*c->getSize();}, dir, 0);
                            vyFluxB += getFlux(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), MOMENTUM, dir, 0);
                            XFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newXx*c->getMass();}, dir, 0);
                            YFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newYx*c->getMass();}, dir, 0);
                            ZFluxB += getFlux2(nextCell, std::abs(vyB*getdt()), children[i]->getWidth(), [](FluidCell *c){return c->newZx*c->getMass();}, dir, 0);
                        }
                        
                    } else {
                        massB = cellB->getMass();
                        EB = cellB->getE(false);
                        if (vyB < 0) {
                            dir = UP;
                            nextCell = cell;
                        } else {
                            dir = DOWN;
                            nextCell = cellB;
                        }
                        massFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c){return c->newMassx;}, dir, 0);
                        EFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c){return c->newEx*c->getSize();}, dir, 0);
                        vyFluxB = getFlux(nextCell, std::abs(vyB*getdt()), cell->getWidth(), MOMENTUM, dir, 0);
                        XFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c){return c->newXx*c->getMass();}, dir, 0);
                        YFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c){return c->newYx*c->getMass();}, dir, 0);
                        ZFluxB = getFlux2(nextCell, std::abs(vyB*getdt()), cell->getWidth(), [](FluidCell *c){return c->newZx*c->getMass();}, dir, 0);
                    }
                    if (vyB < 0) {
                        massFluxB *= -1;
                        EFluxB *= -1;
                        vyFluxB *= -1;
                        XFluxB *= -1;
                        YFluxB *= -1;
                        ZFluxB *= -1;
                    }
                }
                
                cell->newMass = 0.5 * (cell->newMassx + cell->newMassy - massFluxR + massFluxL - massFluxT + massFluxB);
                cell->newE = 0.5 * ((cell->newEx + cell->newEy)*cell->getSize() - EFluxR + EFluxL - EFluxT + EFluxB)/cell->getSize();
            
                double newVx = (mass*cell->getVelocity().getVx() + -vxFluxR + vxFluxL)/cell->newMass;
                double newVy = (mass*cell->getVelocity().getVy() + -vyFluxT + vyFluxB)/cell->newMass;
                
                cell->newVelocity = VelocityVector(newVx, newVy);
                // cell->newX = std::min(1.0f,std::max(0.0f,(float)((mass*X - XFluxR + XFluxL - XFluxT + XFluxB)/cell->newMass)) );
                cell->newX = std::min(1.0f,std::max(0.0f,(float)(0.5 * (cell->newXx*mass + cell->newXy*mass - XFluxR + XFluxL - XFluxT + XFluxB)/cell->newMass)));
                // cell->newY = std::min(1.0f,std::max(0.0f,(float)((mass*Y - YFluxR + YFluxL - YFluxT + YFluxB)/cell->newMass)) );
                cell->newY = std::min(1.0f,std::max(0.0f,(float)(0.5 * (cell->newYx*mass + cell->newYy*mass - YFluxR + YFluxL - YFluxT + YFluxB)/cell->newMass)));
                // cell->newZ = std::min(1.0f,std::max(0.0f,(float)((mass*Z - ZFluxR + ZFluxL - ZFluxT + ZFluxB)/cell->newMass)) );
                cell->newZ = std::min(1.0f,std::max(0.0f,(float)(0.5 * (cell->newZx*mass + cell->newZy*mass - ZFluxR + ZFluxL - ZFluxT + ZFluxB)/cell->newMass)));

            }

            totMass = 0;
            FluidCell *cell;
            for (int i = 0; i < leafCells.size(); i++) {
                cell = leafCells[i];
                if (cell->newMass/cell->getSize() > 100*cell->getDensity()) {
                    // std::cout << cell->depth << " ";
                    // printID(cell->index);
                    // std::cout << std::endl;
                }
                cell->setMass(cell->newMass);
                cell->setE(cell->newE);
                totMass += cell->newMass;
                cell->setVelocity(cell->newVelocity);
                cell->setX(std::max(0.0f,std::min(1.0f,cell->newX)));
                cell->setY(std::max(0.0f,std::min(1.0f,cell->newY)));
                cell->setZ(std::max(0.0f,std::min(1.0f,cell->newZ)));
                cell->sete(cell->getE(false) - 0.5*cell->getDensity()*cell->getVelocity().getMag()*cell->getVelocity().getMag());
                cell->getTemp(true);
                cell->getPressure(true);
            }
        }

        long double getFlux2(FluidCell *cell, long double dist, long double boundary, std::function<long double(FluidCell *)> func, Direction dir, long double totFlux) {
            if (!std::isfinite(dist)) return 0;
            long double leng = (dir == LEFT) || (dir == RIGHT) ? cell->getWidth() : cell->getHeight();
            if (cell->hasChildren()) {
                switch (dir) {
                    case RIGHT:
                        return getFlux2(cell->nw, dist, boundary, func, dir, totFlux) + getFlux2(cell->sw, dist, boundary, func, dir, totFlux);
                        break;
                    case LEFT:
                        return getFlux2(cell->ne, dist, boundary, func, dir, totFlux) + getFlux2(cell->se, dist, boundary, func, dir, totFlux);
                        break;
                    case UP:
                        return getFlux2(cell->se, dist, boundary, func, dir, totFlux) + getFlux2(cell->sw, dist, boundary, func, dir, totFlux);
                        break;
                    case DOWN:
                        return getFlux2(cell->ne, dist, boundary, func, dir, totFlux) + getFlux2(cell->nw, dist, boundary, func, dir, totFlux);
                }
            } else if (dist < leng) {
                return totFlux + func(cell)*boundary*dist/cell->getSize();
            } else {
                FluidCell *nextCell;
                switch (dir) {
                    case RIGHT:
                        nextCell = cell->getRight();
                        break;
                    case LEFT:
                        nextCell = cell->getLeft();
                        break;
                    case UP:
                        nextCell = cell->getTop();
                        break;
                    case DOWN:
                        nextCell = cell->getBottom();
                }
                return getFlux2(nextCell, dist-leng, boundary, func, dir, totFlux + func(cell));
            }
            return 0;
        }

        void setAMR(FluidCell *cell) {
            cell->refinedThisStep = false;
            FluidCell *cellR = cell->getRight();
            FluidCell *cellL = cell->getLeft();
            FluidCell *cellT = cell->getTop();
            FluidCell *cellB = cell->getBottom();
            cellR->absorbChildProps(); cellL->absorbChildProps(); cellT->absorbChildProps(); cellB->absorbChildProps();
            // long double mass = cell->getMass();
            // long double massR = cellR->getMass();
            // long double massL = cellL->getMass();
            // long double massT = cellT->getMass();
            // long double massB = cellB->getMass();
            long double density = cell->getDensity();
            long double densityR = cellR->getDensity();
            long double densityL = cellL->getDensity();
            long double densityT = cellT->getDensity();
            long double densityB = cellB->getDensity();

            if (abs(densityR - density)/density > refineThresh) {
                cell->shouldRefine = true;
                // cellR->shouldRefine = true;
            } else if (abs(densityR - density)/density < coarseThresh) {
                cell->shouldCoarsen = true;
            } else {
                cell->shouldCoarsen = false;
            }
            if (abs(densityL - density)/density > refineThresh) {
                cell->shouldRefine = true;
                // cellL->shouldRefine = true;
            } else if (abs(densityL - density)/density < coarseThresh) {
                cell->shouldCoarsen = true;
            } else {
                cell->shouldCoarsen = false;
            }
            if (abs(densityT - density)/density > refineThresh) {
                // std::cout << "grad " << abs(massT - mass)/mass << std::endl;
                // std::cout << massT << ", " << mass << std::endl;
                cell->shouldRefine = true;
                // cellT->shouldRefine = true;
            } else if (abs(densityT - density)/density < coarseThresh) {
                cell->shouldCoarsen = true;
            } else {
                cell->shouldCoarsen = false;
            }
            if (abs(densityB - density)/density > refineThresh) {
                cell->shouldRefine = true;
                // cellB->shouldRefine = true;
            } else if (abs(densityB - density)/density < coarseThresh) {
                cell->shouldCoarsen = true;
            } else {
                cell->shouldCoarsen = false;
            }

            FluidCell *neighbors[4] = {cellT, cellB, cellL, cellR};
            std::vector<FluidCell *> realNeighbs;
            bool island = true;
            bool willJump = false;
            for (int i = 0; i < 4; i++) {
                if (neighbors[i] != cell) {
                    realNeighbs.push_back(neighbors[i]);
                }
            }
            for (int i = 0; i < realNeighbs.size(); i++) {
                island = island & (int)realNeighbs[i]->hasChildren();
                willJump = willJump | (int)realNeighbs[i]->hasChildren();
            }
            // avoiding flickering of corners
            if (realNeighbs.size() < 4) island = false;
            if (island) {
                cell->shouldRefine = true;
                cell->shouldCoarsen = false;
            } else if (density < densCoarseThresh*maxDensity) {
                cell->shouldRefine = false;
                cell->shouldCoarsen = true;
            } else if (density > densRefineThresh*maxDensity) {
                cell->shouldRefine = true;
                cell->shouldCoarsen = false;
            }
            if (cell->shouldCoarsen) {
                if (cellR->shouldRefine || cellL->shouldRefine || cellT->shouldRefine || cellB->shouldRefine) {
                    cell->shouldCoarsen = false;
                }
                if (willJump) {
                    cell->shouldCoarsen = false;
                }
            }
            if (cell->depth > maxDepth) {
                cell->shouldCoarsen = true;
                cell->shouldRefine = false;
            }

        }

        void checkAMR() {
            std::map<uint64_t,FluidCell*> dictCopy = idToCell;
            FluidCell *cell;
            for (int i = 0; i < leafCells.size(); i++) {
                cell = leafCells[i];
                if (cell->shouldRefine) {
                    TreeLoc loc = getLocFromID(getIDfromLeaf(cell));
                    cell->shouldRefine = false;
                    cell->shouldCoarsen = false;
                    if (loc.depth < maxDepth) {
                        if (!cell->refinedThisStep) refineCell(loc.depth, loc.index);
                    }
                } else if (cell->shouldCoarsen) {
                    uint64_t id = getIDfromLeaf(cell);
                    
                    TreeLoc loc = getLocFromID(id);
                    
                    if (loc.depth > minDepth) {  
                        FluidCell *parent = cell->parent;
                        FluidCell *siblings[4] = {parent->nw, parent->ne, parent->sw, parent->se};
                        bool agreement = true;
                        for (int i = 0; i < 4; i++) {
                            agreement = siblings[i]->shouldCoarsen && agreement;

                            siblings[i]->shouldCoarsen = false;
                            siblings[i]->shouldRefine = false;
                        }
                        
                        if (agreement) {
                            Direction dirs[4] = {DOWN, LEFT, UP, RIGHT};
                            bool resJump = false;
                            for (int i = 0; i < 4; i++) {
                                FluidCell *left = siblings[i]->getLeft();
                                FluidCell *right = siblings[i]->getRight();
                                FluidCell *top = siblings[i]->getTop();
                                FluidCell *bottom = siblings[i]->getBottom();
                                FluidCell *neighbs[4] = {left,right,top,bottom};

                                for (int j = 0; j < 4; j++) {
                                    if (neighbs[j] != cell) {
                                        resJump = resJump | (neighbs[j]->hasChildren() | neighbs[j]->shouldRefine);
                                    }
                                }                         
                            }                            
                            if (!resJump) coarsenCell(loc.depth-1, loc.index >> 2);
                        }
                    }
                    cell->shouldRefine = false;
                    cell->shouldCoarsen = false;
                }
            }
            setNeighbors();
        }

        void updateLeafCells() {
            leafCells.clear();
            for (auto& pair : idToCell) {
                leafCells.push_back(pair.second);
            }
        }

        void update() {
            double new_dt = 0.7*minSize / maxV;
            dt = std::min(new_dt, dt*1.25);
            // std::cout << maxV << ", " << dt << std::endl; 
            // dt = 0.9 * dt + 0.1 * new_dt; 
            maxV = 0;
            auto start = std::chrono::steady_clock::now();
            solveGravPotential(10);
            auto grav = std::chrono::steady_clock::now();
            updateVelocities();
            auto vel = std::chrono::steady_clock::now();
            energyUpdate();
            auto nrg = std::chrono::steady_clock::now();
            advect3();
            auto adv = std::chrono::steady_clock::now();
            // std::cout << "\r" << totMass << std::endl;
            // tote = 0;
            for (int i = 0; i < leafCells.size(); i++) {
                // tote += leafCells[i]->getE(false)*leafCells[i]->getSize();
                setAMR(leafCells[i]);
            }
            // std::cout << "\r" << tote << std::endl;
            auto setamr = std::chrono::steady_clock::now();
            checkAMR();
            auto amr = std::chrono::steady_clock::now();
            updateLeafCells();
            auto leaf = std::chrono::steady_clock::now();

            if (leafCells.size() < 0.5*pow(4,minDepth+2)) maxDepth += 1;
            else if (leafCells.size() > 0.95*pow(4,maxDepth)) maxDepth -= 1;
            // if(std::chrono::duration_cast<std::chrono::milliseconds>(leaf - start).count() > 500) maxDepth -= 1; 

            // std::cout << totMass << std::endl;
            // std::cout << "grav: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(grav - start).count() << "\n";
            // std::cout << "vel: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(vel - grav).count() << "\n";
            // std::cout << "nrg: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(nrg - vel).count() << "\n";
            // std::cout << "adv: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(adv - nrg).count() << "\n";
            // std::cout << "setAMR: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(setamr - adv).count() << "\n";
            // std::cout << "AMR: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(amr - setamr).count() << "\n";
            // std::cout << "leaf: ";
            // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(leaf - amr).count() << "\n";
        }

        void freeGrid() {
            for (int i = 0; i < leafCells.size(); i++) {
                free(leafCells[i]);
            }
        }

        std::vector<FluidCell*> getLeafCells() {
            return leafCells;
        }
        

    private:
        int maxDepth = 7;
        int minDepth = 3;
        // gradient thresh
        float coarseThresh = 0.1; 
        float refineThresh = 0.2;
        // density threshold
        float densCoarseThresh = 0.01;
        float densRefineThresh = 0.5;
        float viscosity = 0;
        double width, height;
        double maxV;
        double minSize; // minimum side length
        float dt;
        FluidCell *root;
        std::map<uint64_t,FluidCell*> idToCell;
        std::vector<FluidCell*> leafCells;
        long double maxDensity;
        long double totMass, tote;
};

class Simulator {
    public:
        Simulator(double width, double height, int startDepth, bool display=true, bool save=false, std::string fn=""): width(width), height(height), display(display), save(save), filename(fn) {
            dt = 0.1;
            // in meters per pixel
            SCALE_H = height / consts::GRID_HEIGHT;
            SCALE_W = width / consts::GRID_WIDTH;
            grid = (FluidGrid *) malloc(sizeof(FluidGrid));
            grid[0] = FluidGrid(width, height, startDepth, dt);
            output.open(filename);
            SDL_Init(SDL_INIT_VIDEO);       // Initializing SDL as Video
            SDL_CreateWindowAndRenderer(consts::GRID_WIDTH, consts::GRID_HEIGHT, 0, &window, &renderer);
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);      // setting draw color
            SDL_RenderClear(renderer);
            
        }

        void drawCells() {
            TreeLoc m;
            uint32_t depth, index;
            xyPos xyBL;
            double cellWidth, cellHeight;
            double thisMaxPressure = 0;
            double thisMaxMass = 0;
            double thisMaxDensity = 0;
            double thisMaxTemperature = 0;
            double thisMinTemperature = 0;
            double thisMinGP = 0;
            double thisMaxe = 0;
            double thisMaxFusion = 0;
            SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255); 
            SDL_RenderClear(renderer);
            std::vector<FluidCell *> leaves = grid->getLeafCells();
            for (int i = 0; i < leaves.size(); i++) {
                
                FluidCell *cell = leaves[i];
                depth = cell->depth;
                index = cell->index;
                double pressure = cell->getPressure(false);
                double density = cell->getDensity();
                double temperature = cell->getTemp(false);
                double gravPotential = cell->getGravPotential();
                double e = cell->getE(false);
                float X = std::max(0.0f,std::min(1.0f,cell->getX()));
                float Y = std::max(0.0f,std::min(1.0f,cell->getY()));
                float Z = std::max(0.0f,std::min(1.0f,cell->getZ()));
                float fusion = cell->getFusionEnergy();
                uint degenerate = cell->getDegenerate();
                xyBL = grid->getXY(depth, index);
                cellWidth = width / pow(pow(4,depth), 0.5);
                cellHeight = height / pow(pow(4,depth), 0.5);

                SDL_Rect rect{(int)(xyBL.x/SCALE_W), (int)(consts::GRID_HEIGHT - ((xyBL.y+cellHeight)/SCALE_H)), (int)(cellWidth/SCALE_W)+1, (int)(cellHeight/SCALE_H)+1};
                if (density > thisMaxDensity) thisMaxDensity = density;
                // if (density > maxDensity) density = maxDensity;
                if (pressure > thisMaxPressure) {
                    thisMaxPressure = pressure;
                }
                if (temperature > thisMaxTemperature ) {
                    // && density > 0.01*maxDensity
                    thisMaxTemperature = temperature;
                }
                if (gravPotential < thisMinGP) thisMinGP = gravPotential;
                if (e > thisMaxe) thisMaxe = e;
                if (fusion > thisMaxFusion) thisMaxFusion = fusion;    

                double maxBit = 255.0;
                double zero = 0;
                if (pressureDisplay) {
                    double scaledP_R = std::min(maxBit,255 * (pressure)/(maxPressure)); 
                    if (scaledP_R < 0) scaledP_R = 0;

                    SDL_SetRenderDrawColor(renderer, scaledP_R, 0, 255, 255);
                } else if (densityDisplay) {                
                    double scaled_dens = std::sqrt(density/maxDensity) * 255;
                    if (scaled_dens > 255) scaled_dens = 255;
                    
                    SDL_SetRenderDrawColor(renderer, scaled_dens*(1-Z), scaled_dens*std::max(0.0f,(1-Y-Z)), scaled_dens*Z, 255);
                    // SDL_RenderFillRect(renderer, &rect);
                } else if (temperatureDisplay) {
                    double scaled_temp = std::min(maxBit,temperature/maxTemperature * 255);
                    // if (density < 0.01*maxDensity) scaled_temp = 0;
                    SDL_SetRenderDrawColor(renderer, 0, scaled_temp, 0, 255);
                } else if (gravPotentialDisplay) {
                    double scaled_pot = std::max(zero,std::min(maxBit,gravPotential/minGP * 255));
                    // if (i == rows/2 && j == cols/2) printf("scaledGP:%f\n",gravPotential/minGP);
                    SDL_SetRenderDrawColor(renderer, 255, scaled_pot, 0, 255);

                } else if (energyDisplay) {
                    double scaled_e = std::max(zero,std::min(maxBit,e/maxe * 255));
                    SDL_SetRenderDrawColor(renderer, scaled_e, 0, 0, 255);
                } else if (hydrogenDisplay) {
                    int scaled_X = std::min(255.0f,X*255);
                    // std::cout << X << std::endl;
                    SDL_SetRenderDrawColor(renderer, scaled_X, scaled_X, 0, 255);
                } else if (heliumDisplay) {
                    int scaled_Y = std::min(255.0f,Y*255);
                    SDL_SetRenderDrawColor(renderer, scaled_Y, 0, scaled_Y, 255);
                } else if (ZDisplay) {
                    int scaled_Z = std::min(255.0f,Z*255);
                    SDL_SetRenderDrawColor(renderer, 0, 0, scaled_Z, 255);   
                } else if (fusionDisplay) {
                    double scaled_fusion = std::max(zero,std::min(maxBit,fusion/maxFusion * 255));
                    if (maxFusion == 0) scaled_fusion = 0;
                    SDL_SetRenderDrawColor(renderer, 0, scaled_fusion, scaled_fusion, 255);
                } else if (degenerateDisplay) {
                    double deg = degenerate*255;
                    SDL_SetRenderDrawColor(renderer, 0, 0, deg, 255);
                }
                
                SDL_RenderFillRect(renderer, &rect);
                if (gridDisplay) {
                    SDL_SetRenderDrawColor(renderer, 255,255,255, 255); // White outline
                    SDL_RenderDrawRect(renderer, &rect);
                }
                
                if (velocityDisplay) {
                    VelocityVector v = cell->getVelocity();
                    Arrow arrow{(int)((xyBL.x+cellWidth/2)/SCALE_W), (int)(consts::GRID_HEIGHT - ((xyBL.y+cellHeight/2)/SCALE_H)), 10*(int)(v.getVx()*grid->getdt()/cellWidth * cellWidth/SCALE_W), -10*(int)(v.getVy()*grid->getdt()/cellHeight * cellHeight/SCALE_H)};
                    SDL_SetRenderDrawColor(renderer, 255,255,255, 255);
                    if (std::isfinite(v.getVx()+v.getVy())) drawArrow(renderer, arrow);
                }
            }
            maxPressure = thisMaxPressure;
            maxDensity = thisMaxDensity;
            // maxDensity = 10;
            maxMass = thisMaxMass;
            maxTemperature = thisMaxTemperature;
            minTemperature = thisMinTemperature;
            // std::cout << maxTemperature << std::endl;
            minGP = thisMinGP; // grav potential
            maxe = thisMaxe;
            maxFusion = thisMaxFusion;
            SDL_RenderPresent(renderer);
            if (thisMaxDensity > maxDensity) maxDensity = thisMaxDensity;
        }
        
        void step(SDL_Event event) {
            grid->update();
            dt = grid->getdt();
            t += grid->getdt();
            if (event.key.type == SDL_KEYDOWN) {
                // printf("key clicked\n");
                switch (event.key.keysym.sym) {
                    case SDLK_p:
                        densityDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        pressureDisplay = !(pressureDisplay);
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        break;
                    case SDLK_d:
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        densityDisplay = !(densityDisplay);
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        std::cout << "max Density: " << maxDensity << std::endl;
                        break;
                    case SDLK_t:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        temperatureDisplay = !(temperatureDisplay);
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        std::cout << "max Temperature: " << maxTemperature << std::endl;
                        break;
                    case SDLK_u:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        energyDisplay = 0;
                        gravPotentialDisplay = !(gravPotentialDisplay);
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        break;
                    case SDLK_e:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = !(energyDisplay);
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        std::cout << "max internal energy: " << maxe << std::endl;
                        break;
                    case SDLK_7:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        hydrogenDisplay = !(hydrogenDisplay);
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        break;
                    case SDLK_8:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        hydrogenDisplay = 0;
                        heliumDisplay = !(heliumDisplay);
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        break;
                    case SDLK_9:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = !(ZDisplay);
                        fusionDisplay = 0;
                        degenerateDisplay = 0;
                        break;
                    case SDLK_f:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = !(fusionDisplay);
                        degenerateDisplay = 0;
                        break;
                    case SDLK_o:
                        densityDisplay = 0;
                        pressureDisplay = 0;
                        temperatureDisplay = 0;
                        gravPotentialDisplay = 0;
                        energyDisplay = 0;
                        hydrogenDisplay = 0;
                        heliumDisplay = 0;
                        ZDisplay = 0;
                        fusionDisplay = 0;
                        degenerateDisplay = !(degenerateDisplay);
                        break;
                    case SDLK_g:
                        gridDisplay = !(gridDisplay);
                        break;
                    case SDLK_v:
                        velocityDisplay = !(velocityDisplay);
                }
            }
            if (display) {
                drawCells();
            }
            if (save) {
                saveSim();
            }
            // SDL_Delay(grid->getdt()*1000);  // setting some Delay
            SDL_Delay(16);
        }

        void saveSim() {
            for (auto &pair : grid->getIDMap()) {
                uint64_t ID = pair.first;
                FluidCell *cell = pair.second;
                long double mass = cell->getMass();
                long double pressure = cell->getPressure(false);
                long double density = cell->getDensity();
                long double temperature = cell->getTemp(false);
                long double gravPotential = cell->getGravPotential();
                long double e = cell->gete(false);
                long double fusion = cell->getFusionEnergy();
                output << "//" << ID << ":" << density << "," << cell->getX() << "," << cell->getY() << "," << cell->getZ() << "," << temperature << "," << pressure << "," << e << "," << gravPotential << "," << fusion;
            }
            output << "\nt:" << t;
            output << "\n";
        }

        void freeSim() {
            grid->freeGrid();
            free(grid);
            output.close();
        }

        double getT() {
            // time passed
            return t;
        }

        FluidGrid *grid;
    private:
        double SCALE_H, SCALE_W;
        double width, height;
        bool display, save;
        double dt;
        SDL_Renderer *renderer = NULL;
        SDL_Window *window = NULL;
        SDL_Surface *screenSurface;
        std::ofstream output = std::ofstream();
        std::string filename = "";
        double t = 0;
        double maxMass = 255;
        double maxPressure = 1;
        double maxDensity = 1;
        double maxTemperature = 1;
        double minTemperature = 0;
        double minGP = 0;
        double maxe = 0;
        double maxFusion = 0;
        uint pressureDisplay = 0; 
        uint temperatureDisplay = 0;
        uint densityDisplay = 1;
        uint gravPotentialDisplay = 0;
        uint energyDisplay = 0;
        uint heliumDisplay = 0;
        uint hydrogenDisplay = 0;
        uint ZDisplay = 0;
        uint fusionDisplay = 0;
        uint degenerateDisplay = 0;
        uint gridDisplay = 0;
        uint velocityDisplay = 0;
};

void printBinaryRecursive(uint64_t num) {
    if (num > 1) {
        printBinaryRecursive(num / 2);
    }
    std::cout << (num % 2);
}


int main(int argv, char **argc) {
    if (argv > 4) std::srand((unsigned) atoi(argc[4]));
    else std::srand((unsigned) std::time(NULL));
    Simulator sim(1e10,1e10,atoi(argc[1]), atoi(argc[2]), atoi(argc[3]), "outputGrid_1e8.txt");
    // Simulator sim(10,10,atoi(argc[1]));

    // for (const auto& pair : sim.grid->getIDMap()) {
    //     printBinaryRecursive(pair.first);
    //     std::cout << ", Value: " << pair.second->getDensity() << std::endl;
    // }
    
    SDL_Event event;
    if (atoi(argc[2])) sim.drawCells();

    SDL_PollEvent(&event);
    while(!(event.type == SDL_QUIT)) {
        SDL_PollEvent(&event);
        // sim.drawCells();
        // if (event.key.state == SDLK_SPACE) {
            sim.step(event);
        // }
        // std::cout << "\r" << sim.getT()/1e8 << std::endl;
        if (atoi(argc[3]) && sim.getT() > 1e8) break;
    }
    sim.freeSim();
    // for (const auto& pair : sim.grid->getIDMap()) {
    //     printBinaryRecursive(pair.first);
    //     std::cout << ", Value: " << pair.second << std::endl;
    // }

    return 0;
}
