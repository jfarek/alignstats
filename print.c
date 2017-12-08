#include "print.h"
#include <string.h>

/**
 * Print percentage of numerator/denomiator to 2 decimal places into buffer.
 */
int print_pct(char *buffer, size_t num, uint64_t numerator, uint64_t denominator)
{
    return snprintf(buffer, num, "%f",
                    (denominator == 0)
                        ? 0.0
                        : (double)numerator / (double)denominator);
                        /*: 100.0 * (double)numerator / (double)denominator);*/
}

/**
 * strncpy() to dest buffer and terminate buffer end with null char.
 */
char *copy_to_buffer(char *dest, const char *src, size_t num)
{
    char *dest_ret;

    dest_ret = strncpy(dest, src, num);
    dest[num - 1] = '\0';

    return dest_ret;
}
