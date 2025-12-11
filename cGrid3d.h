#pragma once 

// a function R3 to R3 represented by Cartesian grid with uniform step
// supports interpolation by means of Laplacian smoothing

class cGrid3d {

public:

struct sPoint {
    Eigen::Vector3d vector;
    bool isPinned = false;
    bool isActive = true;
    Eigen::Vector3d nextVector;
};

void init(const Bounds3d & _bounds, double _step) 
{
    m_step = _step;
    m_bounds = _bounds;

    Eigen::Vector3d span = m_bounds.maxBound - m_bounds.minBound;
    printf("initial span %.6f, %.6f, %.6f\n", span[0], span[1], span[2]);
    
    Ni = static_cast<size_t>(ceil(span[0]/m_step));
    Nj = static_cast<size_t>(ceil(span[1]/m_step)); 
    Nk = static_cast<size_t>(ceil(span[2]/m_step));

    m_vector.resize(Ni*Nj*Nk);
}

void Smooth(int nIter) 
{
    for(int iIter = 0; iIter < nIter; iIter++) {

        // get next position
        for(size_t ii = 0; ii < GetNi(); ii++)
            for(size_t jj = 0; jj < GetNj(); jj++)
                for(size_t kk = 0; kk < GetNk(); kk++) {
                    
                    if(!GetPinned(ii, jj, kk)) {

                        int nAver = 0;
                        Eigen::Vector3d averVec(0.0, 0.0, 0.0);;
                        for(int iC = 0; iC < 3; iC++) {

                            size_t iN = 0, jN = 0, kN = 0;

                            switch(iC) {
                                case 0:
                                    iN = 1;
                                    break;
                                case 1:
                                    jN = 1;
                                    break;
                                case 2:
                                    kN = 2;
                            }

                            if(IsValid(ii+iN, jj+jN, kk+kN) && IsValid(ii-iN, jj-jN, kk-kN)) {
                                averVec += GetVector(ii+iN, jj+jN, kk+kN) + GetVector(ii-iN, jj-jN, kk-kN);
                                nAver += 2;
                            }
                        }

                        if(nAver == 0)
                            GetNextVector(ii, jj, kk) = GetVector(ii, jj, kk);
                        else
                            GetNextVector(ii, jj, kk) = averVec * (1.0/nAver);
                    }
                }

        // set next position        
        for(size_t ii = 0; ii < GetNi(); ii++)
            for(size_t jj = 0; jj < GetNj(); jj++)
                for(size_t kk = 0; kk < GetNk(); kk++) {
                        
                    if(!GetPinned(ii, jj, kk)) {
                        GetVector(ii, jj, kk) = GetNextVector(ii, jj, kk);
                    }
                }
    }
}


// size_t GetN() const { return N; }
size_t GetNi() const { return Ni; }
size_t GetNj() const { return Nj; }
size_t GetNk() const { return Nk; }

Eigen::Vector3d GetPosition(size_t i, size_t j, size_t k) { 
    Eigen::Vector3d shift = m_bounds.minBound; 
    shift[0] += i*m_step; shift[1] += j*m_step; shift[2] += k*m_step;
    return shift;
}

Eigen::Vector3d& GetVector(size_t i) { return m_vector[i].vector; }
Eigen::Vector3d& GetVector(size_t i, size_t j, size_t k) { return m_vector[getIndex(i, j, k)].vector; }

Eigen::Vector3d& GetNextVector(size_t i) { return m_vector[i].nextVector; }
Eigen::Vector3d& GetNextVector(size_t i, size_t j, size_t k) { return m_vector[getIndex(i, j, k)].nextVector; }

bool& GetPinned(size_t i) { return m_vector[i].isPinned; }
bool& GetPinned(size_t i, size_t j, size_t k) { return m_vector[getIndex(i, j, k)].isPinned; }

bool IsValid(size_t i, size_t j, size_t k) const { return i < Ni && j < Nj && k < Nk; }
size_t getIndex(size_t i, size_t j, size_t k) const { return i + Ni*j + Ni*Nj*k;}
double getStep() const { return m_step;}
const Bounds3d & getBounds() const {return m_bounds;}

private:

    size_t Ni = 0;
    size_t Nj = 0;
    size_t Nk = 0;

    std::vector<sPoint> m_vector;

    double   m_step;
    Bounds3d m_bounds;
};

