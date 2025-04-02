/*-----------------------------------------------------------------------------
 * Planck
 * Copyright (C) 2024 Hemanth Haridas, University of Utah
 * Contact: hemanthhari23@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or a later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 ----------------------------------------------------------------------------*/

#include "planck_formdump.h"

void formatBinaryDump(std::fstream *filePointer, std::fstream *outPointer, std::error_code *errorFlag, std::string *errorMessage)
{
    // first reset the error buffers
    errorFlag->clear();
    errorMessage->clear();

    // check if the filepointer is valid
    if (!filePointer || !filePointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable To Open The Input File. Please Check The Input File Provided.";
        return;
    }

    // Determine the size of the file
    filePointer->seekg(0, std::ios::end); // Move to the end of the file
    std::streampos fileSize = filePointer->tellg(); // Get the file size
    filePointer->seekg(0, std::ios::beg); // Move back to the start of the file

    std::vector<std::double_t> fileContents(fileSize / sizeof(std::double_t));
    filePointer->read(reinterpret_cast<char*>(fileContents.data()), fileSize);

    // check if the filepointer is valid
    if (!outPointer || !outPointer->is_open())
    {
        *errorFlag = std::make_error_code(std::errc::io_error);
        *errorMessage = "Unable To Open The Formatted File For Writing.";
        return;
    }

    // Now write the header section
    // std::ostringstream oss;
    // oss << "Number of Basis Functions : " << static_cast<int>(sqrt(fileContents.size()));
    // std::string headerLine = oss.str();
    // std::cout << headerLine.size() << "\n";
    // outPointer->write(headerLine.c_str(), headerLine.size());

    // Now write the matrix
    // outPointer->close();
}

int main(int argc, char const *argv[])
{
    std::error_code error_flag;
    std::string error_message;

    // check if an input file is provided. Exit if not input file is provided
    if (argc < 2)
    {
        std::cout << std::setw(20) << std::left << "[Error]   <= " << std::left << " Unable To Find An Input File. Please Run formDump As : formDump input " << "\n";
        exit(-1);
    }

    // proceed to parsing to the input file
    std::string input_file = argv[1];
    std::string output_file = input_file + ".formdump";
    std::fstream file_pointer(input_file);
    std::fstream out_pointer(output_file, std::ios::out | std::ios::trunc);

    formatBinaryDump(&file_pointer, &out_pointer, &error_flag, &error_message);

    // check if input file was parsed correctly
    if (error_flag && error_flag.value() != 71)
    {
        std::cout << std::setw(21) << std::left << "[Error]   <=  " << std::left << error_message << "\n";
        exit(error_flag.value());
    }

    exit(0);
}