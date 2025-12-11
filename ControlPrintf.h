#pragma once

#include <iostream>
#include <string>
#include <sstream>

// define a prinf named instance of this class to disable printf's 

class ControlPrintf {
public:

    ControlPrintf(bool isActive = false) : _isActive(isActive) {}

    // Define the operator() as a template function
    template <typename... Args>
    void operator()(const char* str, Args&&... args) const {
	    if(_isActive)
            printf(str, std::forward<Args>(args)...);
    }

private:
    bool _isActive = false;
};

