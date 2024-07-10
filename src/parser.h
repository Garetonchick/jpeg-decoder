#pragma once
#include <iostream>
#include <vector>
#include "stream_reader.h"
#pragma once
#include "jpeg.h"

JPEGData Parse(std::istream& in);
std::string ParseResultLogString(const JPEGData& data);