#!/usr/bin/awk -f

BEGIN {
    RS = "\n"
    file_counter = 1
    output_file = "6-31g-" file_counter ".basis"
}

{
    if ($0 == "****") {
        close(output_file)
        file_counter++
        output_file = "6-31g-" file_counter ".basis"
    } else {
        print $0 >> output_file
    }
}