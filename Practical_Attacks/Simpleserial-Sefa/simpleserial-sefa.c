/*
    This file is part of the ChipWhisperer Example Targets
    Copyright (C) 2012-2017 NewAE Technology Inc.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdlib.h>

uint8_t get_mask(uint8_t *m, uint8_t len)
{
    return 0x00;
}

uint8_t get_key(uint8_t *k, uint8_t len)
{
    return 0x00;
}

uint8_t get_pt(uint8_t *pt, uint8_t len)
{
    uint8_t plane1[5];
    uint8_t plane2[5];
    uint8_t randomness[5];

    for (size_t i = 0; i < 5; i++)
    {
        plane1[i] = pt[i + 0];
        plane2[i] = pt[i + 5];
        randomness[i] = pt[i + 10];
    }

    // ////////////////////////////////
    // //      MASKED CHI - DOM      //
    // ////////////////////////////////

    trigger_high();

    asm volatile(
        "mov r20, %[a0]"
        "\n\t"
        "com r20"
        "\n\t"
        "mov r21, %[b0]"
        "\n\t"
        "and r21, r20"
        "\n\t"
        "mov r22, %[b1]"
        "\n\t"
        "and r22, r20"
        "\n\t"
        "mov r23, %[a1]"
        "\n\t"
        "and r23, %[b0]"
        "\n\t"
        "mov r24, %[b1]"
        "\n\t"
        "and r24, %[a1]"
        "\n\t"
        "eor r22, %[R0]"
        "\n\t"
        "eor r23, %[R0]"
        "\n\t"
        "eor r21, r22"
        "\n\t"
        "eor r23, r24"
        "\n\t"
        "mov r25, %[e0]"
        "\n\t"
        "mov r26, %[e1]"
        "\n\t"
        "eor r25, r21"
        "\n\t"
        "eor r26, r23"
        "\n\t"

        //---------------------

        "mov r20, %[b0]"
        "\n\t"
        "com r20"
        "\n\t"
        "mov r21, %[c0]"
        "\n\t"
        "and r21, r20"
        "\n\t"
        "mov r22, %[c1]"
        "\n\t"
        "and r22, r20"
        "\n\t"
        "mov r23, %[b1]"
        "\n\t"
        "and r23, %[c0]"
        "\n\t"
        "mov r24, %[c1]"
        "\n\t"
        "and r24, %[b1]"
        "\n\t"
        "eor r22, %[R1]"
        "\n\t"
        "eor r23, %[R1]"
        "\n\t"
        "eor r21, r22"
        "\n\t"
        "eor r23, r24"
        "\n\t"
        "mov r18, %[a0]"
        "\n\t"
        "mov r19, %[a1]"
        "\n\t"
        "eor r18, r21"
        "\n\t"
        "eor r19, r23"
        "\n\t"

        //---------------------

        "mov r20, %[c0]"
        "\n\t"
        "com r20"
        "\n\t"
        "mov r21, %[d0]"
        "\n\t"
        "and r21, r20"
        "\n\t"
        "mov r22, %[d1]"
        "\n\t"
        "and r22, r20"
        "\n\t"
        "mov r23, %[c1]"
        "\n\t"
        "and r23, %[d0]"
        "\n\t"
        "mov r24, %[d1]"
        "\n\t"
        "and r24, %[c1]"
        "\n\t"
        "eor r22, %[R2]"
        "\n\t"
        "eor r23, %[R2]"
        "\n\t"
        "eor r21, r22"
        "\n\t"
        "eor r23, r24"
        "\n\t"
        "eor %[b0], r21"
        "\n\t"
        "eor %[b1], r23"
        "\n\t"

        //---------------------

        "mov r20, %[d0]"
        "\n\t"
        "com r20"
        "\n\t"
        "mov r21, %[e0]"
        "\n\t"
        "and r21, r20"
        "\n\t"
        "mov r22, %[e1]"
        "\n\t"
        "and r22, r20"
        "\n\t"
        "mov r23, %[d1]"
        "\n\t"
        "and r23, %[e0]"
        "\n\t"
        "mov r24, %[e1]"
        "\n\t"
        "and r24, %[d1]"
        "\n\t"
        "eor r22, %[R3]"
        "\n\t"
        "eor r23, %[R3]"
        "\n\t"
        "eor r21, r22"
        "\n\t"
        "eor r23, r24"
        "\n\t"
        "eor %[c0], r21"
        "\n\t"
        "eor %[c1], r23"
        "\n\t"

        //---------------------

        "com %[e0]"
        "\n\t"
        "mov r21, %[a0]"
        "\n\t"
        "and r21, %[e0]"
        "\n\t"
        "mov r22, %[a1]"
        "\n\t"
        "and r22, %[e0]"
        "\n\t"
        "mov r23, %[e1]"
        "\n\t"
        "and r23, %[a0]"
        "\n\t"
        "mov r24, %[a1]"
        "\n\t"
        "and r24, %[e1]"
        "\n\t"
        "eor r22, %[R4]"
        "\n\t"
        "eor r23, %[R4]"
        "\n\t"
        "eor r21, r22"
        "\n\t"
        "eor r23, r24"
        "\n\t"
        "eor %[d0], r21"
        "\n\t"
        "eor %[d1], r23"
        "\n\t"

        //---------------------

        "mov %[a0], r18"
        "\n\t"
        "mov %[a1], r19"
        "\n\t"
        "mov %[e0], r25"
        "\n\t"
        "mov %[e1], r26"
        "\n\t"
        : "=r"(plane1[0]), "=r"(plane2[0]),
          "=r"(plane1[1]), "=r"(plane2[1]),
          "=r"(plane1[2]), "=r"(plane2[2]),
          "=r"(plane1[3]), "=r"(plane2[3]),
          "=r"(plane1[4]), "=r"(plane2[4])
        :
        [a0] "0"(plane1[0]), [a1] "1"(plane2[0]),
        [b0] "2"(plane1[1]), [b1] "3"(plane2[1]),
        [c0] "4"(plane1[2]), [c1] "5"(plane2[2]),
        [d0] "6"(plane1[3]), [d1] "7"(plane2[3]),
        [e0] "8"(plane1[4]), [e1] "9"(plane2[4]),
        [R0] "r"(randomness[0]),
        [R1] "r"(randomness[1]),
        [R2] "r"(randomness[2]),
        [R3] "r"(randomness[3]),
        [R4] "r"(randomness[4])
        : "r20", "r21", "r22", "r23", "r24", "r25", "r26", "r18", "r19");

    // ////////////////////////////////////
    // //      MASKED CHI - TOFFOLI      //
    // ////////////////////////////////////
    // uint8_t r0 = randomness[0];
    // uint8_t r1 = r0;

    // trigger_high();

    // asm volatile(
    //     "mov r24, %[e0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "mov r25, %[a1]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[r0], r25"
    //     "\n\t"

    //     "mov r24, %[e0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[a0]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[r0], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[a1]"
    //     "\n\t"
    //     "and r25, %[e1]"
    //     "\n\t"
    //     "eor %[r1], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[a0]"
    //     "\n\t"
    //     "and r25, %[e1]"
    //     "\n\t"
    //     "eor %[r1], r25"
    //     "\n\t"

    //     //---------------------

    //     "mov r24, %[b0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "mov r25, %[c1]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[a0], r25"
    //     "\n\t"

    //     "mov r24, %[b0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[c0]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[a0], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[c1]"
    //     "\n\t"
    //     "and r25, %[b1]"
    //     "\n\t"
    //     "eor %[a1], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[c0]"
    //     "\n\t"
    //     "and r25, %[b1]"
    //     "\n\t"
    //     "eor %[a1], r25"
    //     "\n\t"

    //     //---------------------

    //     "mov r24, %[d0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "mov r25, %[e1]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[c0], r25"
    //     "\n\t"

    //     "mov r24, %[d0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[e0]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[c0], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[e1]"
    //     "\n\t"
    //     "and r25, %[d1]"
    //     "\n\t"
    //     "eor %[c1], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[e0]"
    //     "\n\t"
    //     "and r25, %[d1]"
    //     "\n\t"
    //     "eor %[c1], r25"
    //     "\n\t"

    //     //---------------------

    //     "mov r24, %[a0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "mov r25, %[b1]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[e0], r25"
    //     "\n\t"

    //     "mov r24, %[a0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[b0]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[e0], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[b1]"
    //     "\n\t"
    //     "and r25, %[a1]"
    //     "\n\t"
    //     "eor %[e1], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[b0]"
    //     "\n\t"
    //     "and r25, %[a1]"
    //     "\n\t"
    //     "eor %[e1], r25"
    //     "\n\t"

    //     //---------------------

    //     "mov r24, %[c0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "mov r25, %[d1]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[b0], r25"
    //     "\n\t"

    //     "mov r24, %[c0]"
    //     "\n\t"
    //     "com r24"
    //     "\n\t"
    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[d0]"
    //     "\n\t"
    //     "and r25,   r24"
    //     "\n\t"
    //     "eor %[b0], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[d1]"
    //     "\n\t"
    //     "and r25, %[c1]"
    //     "\n\t"
    //     "eor %[b1], r25"
    //     "\n\t"

    //     "ldi r25, 0x00"
    //     "\n\t"
    //     "mov r25, %[d0]"
    //     "\n\t"
    //     "and r25, %[c1]"
    //     "\n\t"
    //     "eor %[b1], r25"
    //     "\n\t"

    //     //---------------------

    //     "eor %[d0], %[r0]"
    //     "\n\t"
    //     "eor %[d1], %[r1]"
    //     "\n\t"

    //     //---------------------

    //     : "=r"(plane1[0]), "=r"(plane2[0]),
    //       "=r"(plane1[1]), "=r"(plane2[1]),
    //       "=r"(plane1[2]), "=r"(plane2[2]),
    //       "=r"(plane1[3]), "=r"(plane2[3]),
    //       "=r"(plane1[4]), "=r"(plane2[4]),
    //       "=r"(r0), "=r"(r1)
    //     :
    //     [a0] "0"(plane1[0]), [a1] "1"(plane2[0]),
    //     [b0] "2"(plane1[1]), [b1] "3"(plane2[1]),
    //     [c0] "4"(plane1[2]), [c1] "5"(plane2[2]),
    //     [d0] "6"(plane1[3]), [d1] "7"(plane2[3]),
    //     [e0] "8"(plane1[4]), [e1] "9"(plane2[4]),
    //     [r0] "10"(r0), [r1] "11"(r1)
    //     : "r24", "r25" // t0, t1
    // );

    trigger_low();

    for (size_t i = 0; i < 5; i++)
    {
        pt[i + 0] = plane1[i];
        pt[i + 5] = plane2[i];
        pt[i + 10] = 0;
    }

    simpleserial_put('r', 16, pt);
    return 0x00;
}

uint8_t reset(uint8_t *x, uint8_t len)
{
    // Reset key here if needed
    return 0x00;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();

    /*
    putch('h');
    putch('e');
    putch('l');
    putch('l');
    putch('o');
    putch('\n');
    */

    simpleserial_init();
    simpleserial_addcmd('k', 16, get_key);
    simpleserial_addcmd('p', 16, get_pt);
    simpleserial_addcmd('x', 0, reset);
    while (1)
        simpleserial_get();
}
