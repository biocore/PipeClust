#ifndef __UTIL_H__
#define __UTIL_H__

// Error levels
#define INFO_MSG 0
#define WARN_ERROR 1
#define FATAL_ERROR 2

/*
    Prints an error through stderr
    

    Inputs:
        level: integer that indicates if the level of the error
        format: string that contains the text to be written to the stream.
            It can optionally contain embedded format specifiers that are
            replaced by the values specified in subsequent additional
            arguments and formatted as requested.
*/
void error_handler(int level, char* format, ...);

#endif