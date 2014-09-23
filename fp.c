#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

/*
    Bit Fields: S = sign, E = exponent, M = Mantissa
    S E E E  E E E E  E M M M  M M M M  M M M M  M M M M  M M M M  M M M M
   --|-----------------|---------------------------------------------------
    1 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0
    8        0        0        0        0        0        0        0        = Sign Bitmask

    0 1 1 1  1 1 1 1  1 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0 
    7        F        8        0        0        0        0        0        = Exponent Bitmask

    0 0 0 0  0 0 0 0  0 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1  1 1 1 1
    0        0        7        F        F        F        F        F        = Mantissa Bitmask
*/

#define FLIB_SIGN_MASK      (0x80000000)
#define FLIB_EXPONENT_MASK  (0x7F800000)
#define FLIB_MANTISSA_MASK  (0x007FFFFF)
#define FLIB_CARRY_MASK     (0xFF800000)

#define FLIB_SIGN_SHIFT     (31)
#define FLIB_EXPONENT_SHIFT (23)
#define FLIB_MANTISSA_SHIFT (0)
#define FLIB_CARRY_SHIFT    (23)

#define FLIB_MASKOFF(v,m)        ((v) & (m))
#define FLIB_EXTRACT_PART(v,m,s) (FLIB_MASKOFF((v),(m)) >> (s))
#define FLIB_PUSH_PART(v,m,s)    (FLIB_MASKOFF(((v)<<(s)), (m)))

typedef struct {
    uint8_t sign;
    int16_t exponent;
    uint32_t mantissa;
    uint32_t characteristic;
} flib_float_t;


void flib_explode(float f, flib_float_t* s);
float flib_pack(flib_float_t* s);
void flib_adjust_radix(flib_float_t* s, int v, bool round);
void flib_align_radix(flib_float_t* ft_a, flib_float_t* ft_b);
void flib_correct_radix(flib_float_t* ft_r);
float __aeabi_fsub(float a, float b);
float __aeabi_fsub_struct(flib_float_t* ft_a, flib_float_t* ft_b);
float __aeabi_fadd(float a, float b);
float __aeabi_fadd_struct(flib_float_t* ft_a, flib_float_t* ft_b);
float __aeabi_fmul(float a, float b);
float __aeabi_fmul_struct(flib_float_t* ft_a, flib_float_t* ft_b);
float __aeabi_fdiv(float a, float b);
float __aeabi_fdiv_struct(flib_float_t* ft_a, flib_float_t* ft_b);
float __aeabi_ui2f(unsigned int i);

void ftoa(float f, char* buffer);
float atof(char* buffer);

/** Assist function. Problem case was no dependency on built-ins **/
uint32_t __uipow(uint32_t  base, uint32_t exponent) {
    if(exponent == 0) return 1;
    return base * __uipow(base, --exponent);
}

void print_float_t(flib_float_t* s) {
    //DEBUGGING
    printf("%s%u.",(s->sign == 0 ? "" : "-"), s->characteristic);
    uint32_t m, bit;
    for(m = FLIB_MASKOFF(s->mantissa, FLIB_MANTISSA_MASK); FLIB_MASKOFF(m, FLIB_MANTISSA_MASK) > 0; m <<= 1) {
        if((m & 0x400000) > 0)
            putchar('1');
        else
            putchar('0');
    }
    printf(" x 2^(%i)\n", s->exponent - 127);
}

void print_float(float f) {
    flib_float_t ft;
    flib_explode(f, &ft);
    print_float_t(&ft);
}


void flib_explode(float f, flib_float_t* s) {
    uint32_t *f_raw = (uint32_t*)&f;
    s->sign     = FLIB_EXTRACT_PART(*f_raw, FLIB_SIGN_MASK, FLIB_SIGN_SHIFT);
    s->exponent = FLIB_EXTRACT_PART(*f_raw, FLIB_EXPONENT_MASK, FLIB_EXPONENT_SHIFT);
    s->mantissa = FLIB_EXTRACT_PART(*f_raw, FLIB_MANTISSA_MASK, FLIB_MANTISSA_SHIFT);
    s->characteristic = (*f_raw > 0) ? 1 : 0;
}

float flib_pack(flib_float_t* s) {
    float _ret;
    uint32_t *f_raw = (uint32_t*)&_ret;
    flib_correct_radix(s);
    *f_raw = FLIB_PUSH_PART(s->sign, FLIB_SIGN_MASK, FLIB_SIGN_SHIFT) |
             FLIB_PUSH_PART(s->exponent, FLIB_EXPONENT_MASK, FLIB_EXPONENT_SHIFT) |
             FLIB_PUSH_PART(s->mantissa, FLIB_MANTISSA_MASK, FLIB_MANTISSA_SHIFT);
    return _ret;
}

void flib_adjust_radix(flib_float_t* s, int v, bool round) {
    
    /* The variable v represents the direction to move the radix point. This, a value of 2 will     *
     * move the radix 2 places to the left and a value of -2 will move the radix right.             *
     *                                                                                              *
     * One edge case is a round operation. This only occurse when moving the radix to the left.     *
     * In this event it is important to only move left if requiered and exectute a left             *
     * move only 1 time per operation. This code will store the dropped bits and check              *
     * if a round is needed. It is possible that by rounding up you could cause a carry             *
     * from the mantissa to the characteristic. If this happens you MUST call this function         *
     * again. Rounding should not happen the second round since the carry would have forced a 0     *
     * into the last bit.                                                                           */

    #define __FLIB_LAST_BIT  (0x01)
    #define __FLIB_TOP_BIT   (0x00400000)
    #define __FLIB_ROUND_BIT (0x80000000)
    #define __FLIB_MANTISSA_TOP_SHIFT (22)

    uint32_t tmp, round_bits = 0;
    s->exponent += v;
    if(v > 0) {
        for(;v>0;v--) {

            // 1.01010 x 2^2
            // 0.10101 x 2^3

            if(round && FLIB_MASKOFF(s->mantissa, __FLIB_LAST_BIT) == 1)
                round_bits = ((round_bits >> 1) | __FLIB_ROUND_BIT);

            else if(round)
                round_bits = (round_bits >> 1);

            tmp = FLIB_MASKOFF(s->characteristic, __FLIB_LAST_BIT);
            s->characteristic >>= 1;
            s->mantissa = (s->mantissa >> 1) | (tmp << __FLIB_MANTISSA_TOP_SHIFT);
        }

        if(round && FLIB_MASKOFF(round_bits, __FLIB_ROUND_BIT)) {
            s->mantissa += 1;
            s->characteristic += FLIB_EXTRACT_PART(s->mantissa, FLIB_CARRY_MASK, FLIB_CARRY_SHIFT);
        }

    } else if(v < 0) {
        for(;v<0;v++) {

            // 1.01010 x 2^2
            // 10.1010 x 2^1

            tmp = FLIB_MASKOFF(s->mantissa, __FLIB_TOP_BIT) >> __FLIB_MANTISSA_TOP_SHIFT;
            s->characteristic = (s->characteristic << 1) | (tmp & __FLIB_LAST_BIT);
            s->mantissa = FLIB_MASKOFF((s->mantissa << 1), FLIB_MANTISSA_MASK);
        }
    }
}

void flib_align_radix(flib_float_t* ft_a, flib_float_t* ft_b) {

   /* We prefer to only shift in a non distructive direction       *
    * But we can only shift that way 31 times. We have to catch    *
    * any shifts more then 31 and shift the remaining amount       *
    * in a distructive direction on the other addend/augend        */

    #define __FLIB_ADJUST_LIMIT (32)
    #define __FLIB_ADJUST_MAX   (31)

    int diff = ft_a->exponent - ft_b->exponent;
    if(diff > 0) {

        if(diff < __FLIB_ADJUST_LIMIT)
            // Safe direction, no carry
            flib_adjust_radix(ft_a, -diff, false);

        else {
            // Safe for a but not for b, carry b
            flib_adjust_radix(ft_a, -(__FLIB_ADJUST_MAX), false);
            flib_adjust_radix(ft_b, diff-__FLIB_ADJUST_MAX, true);
        }

    }
    else if(diff < 0) {

        if(diff > -__FLIB_ADJUST_LIMIT)
            // Safe direction, no carry
            flib_adjust_radix(ft_b, diff, false);

        else {
            // Safe for b but not for a, carry a
            flib_adjust_radix(ft_b, -(__FLIB_ADJUST_MAX), false);
            flib_adjust_radix(ft_a, -(diff+__FLIB_ADJUST_MAX), true);
        }

    }
}

void flib_correct_radix(flib_float_t* ft_r) {

    /* If the radix is moving in a safe direction, ie: a non distructive direction,     *
     * we can simply loop untill it is correct.                                         *
     * If the radix is moving in a unsafe direction,                                    *
     * We must calculate the move count to properly round the bits that fall off.       *
     * An additional test it done at the end to see if rounding caused a larger carry.  *
     *                                                                                  *
     * One edge case is a zero value. Capture this one explicitly and set the exponent  */
     
    uint32_t tmp, f_adjust;
    if(ft_r->characteristic > 1) {

        // to round properly, we need to do move the radix with 1 call to fadjust.
        // Count how many bits we need to shift. 
        // This for loop is f(x) = ceil(log(x)/log(2)) - 1
        for(tmp = ft_r->characteristic, f_adjust = -1; tmp > 0; tmp>>=1, f_adjust++);
        flib_adjust_radix(ft_r, f_adjust, true);

        // It is possible that a round could cause us to need 1 more shift. 
        // If this is the case, rounding will make no diffrence. 
        // Example: 1111.1111 to shift left by 3 produces 10.0. 
        //          10.0 will not be a rounding issue.
        if(ft_r->characteristic > 1) 
            flib_adjust_radix(ft_r, 1, false);

    } else if (ft_r->characteristic == 0 && ft_r->mantissa == 0) {

        // Edge case
        ft_r->exponent = 0;
        ft_r->sign = 0; // can't have negative zero

    } else{

        // Shifting in a safe direction. just roll through.
        while(ft_r->characteristic < 1 && ft_r->mantissa != 0)
            flib_adjust_radix(ft_r, -1, false);

    }
}

float __aeabi_fsub(float a, float b) {
    flib_float_t ft_a;
    flib_float_t ft_b;
    flib_explode(a, &ft_a);
    flib_explode(b, &ft_b);
    return __aeabi_fsub_struct(&ft_a, &ft_b);
}

float __aeabi_fsub_struct(flib_float_t* ft_a, flib_float_t* ft_b) {
    
    #define __FLIB_FORCED_CARRY_FLAG (0x00800000)

    flib_float_t ft_r, *ft_ra, *ft_rb;
    uint8_t sign;

    if(ft_a->sign != ft_b->sign && ft_a->sign == 1) {
        ft_b->sign = 1;
        return __aeabi_fadd_struct(ft_a, ft_b);
    } else if(ft_a->sign != ft_b->sign) {
        ft_b->sign = 0;
        return __aeabi_fadd_struct(ft_a, ft_b);
    }

    flib_align_radix(ft_a, ft_b);

    // Always subtract FROM the LARGER value. 
    // If we have to swap the subtrahend and the minuend then the sign must be flipped.
    if( (ft_b->characteristic > ft_a->characteristic) ||
        (ft_b->characteristic == ft_a->characteristic && ft_b->mantissa > ft_a->mantissa) ) {
        // Sign Flip
        sign = ft_a->sign > 0 ? 0 : 1;
        ft_ra = ft_b;
        ft_rb = ft_a;
    } else {
        sign = ft_a->sign > 0 ? 1 : 0;
        ft_ra = ft_a;
        ft_rb = ft_b;
    }

    // Check for barrow
    if(ft_rb->mantissa > ft_ra->mantissa) {
        ft_ra->characteristic--;
        ft_ra->mantissa |= __FLIB_FORCED_CARRY_FLAG; // Stick in the cary bit.
    }

    ft_r.characteristic = ft_ra->characteristic - ft_rb->characteristic;
    ft_r.mantissa = FLIB_MASKOFF((ft_ra->mantissa - ft_rb->mantissa), FLIB_MANTISSA_MASK);
    ft_r.sign = sign;
    ft_r.exponent = ft_ra->exponent;

    return flib_pack(&ft_r);
}

float __aeabi_fadd(float a, float b) {
    flib_float_t ft_a;
    flib_float_t ft_b;
    flib_explode(a, &ft_a);
    flib_explode(b, &ft_b);
    return __aeabi_fadd_struct(&ft_a, &ft_b);
}

float __aeabi_fadd_struct(flib_float_t* ft_a, flib_float_t* ft_b) {
    flib_float_t ft_r;

    if(ft_a->sign != ft_b->sign && ft_a->sign == 1) {
        ft_a->sign = 0;
        return __aeabi_fsub_struct(ft_b, ft_a);
    } else if(ft_a->sign != ft_b->sign) {
        ft_b->sign = 0;
        return __aeabi_fsub_struct(ft_a, ft_b);
    }

    flib_align_radix(ft_a, ft_b);

    uint32_t n_mantissa = ft_a->mantissa + ft_b->mantissa;
    uint32_t a_carry = FLIB_EXTRACT_PART(n_mantissa, FLIB_CARRY_MASK, FLIB_CARRY_SHIFT);

    ft_r.mantissa = FLIB_MASKOFF(n_mantissa, FLIB_MANTISSA_MASK);
    ft_r.characteristic = ft_a->characteristic + ft_b->characteristic + a_carry;
    ft_r.sign = ft_a->sign;
    ft_r.exponent = ft_a->exponent;

    return flib_pack(&ft_r);
}

float __aeabi_fmul(float a, float b) {
    flib_float_t ft_a;
    flib_float_t ft_b;
    flib_explode(a, &ft_a);
    flib_explode(b, &ft_b);
    return __aeabi_fmul_struct(&ft_a, &ft_b);
}

float __aeabi_fmul_struct(flib_float_t* ft_a, flib_float_t* ft_b) {
    flib_float_t ft_r;
    
    int n_exponent = (
        (ft_a->exponent - 127) + 
        (ft_b->exponent - 127)
    ) + 127;

    int n_mantissa = (int)((
        (long long unsigned int)(ft_a->mantissa | (ft_a->characteristic << 23)) * 
        (long long unsigned int)(ft_b->mantissa |  (ft_a->characteristic << 23))
    ) >> 23);

    ft_r.characteristic = FLIB_EXTRACT_PART(n_mantissa, FLIB_CARRY_MASK, FLIB_CARRY_SHIFT);
    ft_r.mantissa = FLIB_MASKOFF(n_mantissa, FLIB_MANTISSA_MASK);
    ft_r.exponent = n_exponent;
    ft_r.sign = (ft_a->sign + ft_b->sign) % 2;

    return flib_pack(&ft_r);
}

float __aeabi_fdiv(float a, float b) {
    flib_float_t ft_a;
    flib_float_t ft_b;
    flib_explode(a, &ft_a);
    flib_explode(b, &ft_b);
    return __aeabi_fdiv_struct(&ft_a, &ft_b);
}

/* This function uses reciprocal division. It is much faster but sligtly less accurate.
 * This method also requires Uldiv support. To remove dependencies i have replaced
 * it with the function below that uses long division.
 */ /*
float __aeabi_fdiv_struct(flib_float_t* ft_a, flib_float_t* ft_b) {
    flib_float_t ft_r;

    // using reciprical division

    if (ft_b->mantissa == 0 && ft_b->characteristic == 0 && ft_b->exponent == 0) {

        // This is infinity. It is a special value in FP structure.

        ft_r.characteristic = 1;
        ft_r.mantissa = 0;
        ft_r.exponent = 255;
        return flib_pack(&ft_r);
    }

    uint64_t dividend = 1ULL<<46;
    uint64_t divisor = (uint64_t)ft_b->mantissa;

    divisor |= (ft_b->characteristic << 23);
    uint32_t quotient = (uint32_t)(dividend / divisor);

    ft_r.characteristic = FLIB_EXTRACT_PART(quotient, FLIB_CARRY_MASK, FLIB_CARRY_SHIFT);
    ft_r.mantissa = FLIB_MASKOFF(quotient, FLIB_MANTISSA_MASK);
    ft_r.exponent = (-(ft_b->exponent - 127)) + 127;
    ft_r.sign = ft_b->sign;

    // Multiply expects a correct form so fix the radix

    flib_correct_radix(&ft_r);
    return __aeabi_fmul_struct(ft_a, &ft_r);
}
*/

float __aeabi_fdiv_struct(flib_float_t* ft_a, flib_float_t* ft_b) {
    flib_float_t ft_r;

    uint64_t dividend = ft_a->mantissa | (ft_a->characteristic << 23);
    uint64_t divisor = ft_b->mantissa | (ft_b->characteristic << 23);
    uint32_t quotient = 0;

    int n_exponent = (
        (ft_a->exponent - 127) -
        (ft_b->exponent - 127)
    ) + 127;

    if (ft_b->mantissa == 0 && ft_b->characteristic == 0 && ft_b->exponent == 0) {
        /** This is infinity. It is a special value in FP structure. **/
        ft_r.characteristic = 1;
        ft_r.mantissa = 0;
        ft_r.exponent = 255;
        return flib_pack(&ft_r);
    }

    // Increase the percision but make sure divisor is greator then dividend
    dividend <<= 31;
    divisor <<= 32;

    // This will messup my exponent. 
    // So wee need to track them to subtract off the exponent after the divide.
    int skips = 0;
    int skip_flag = 0;

    for(;divisor > 0 && quotient < (FLIB_MANTISSA_MASK | 0x800000); divisor>>=1) {
        if(divisor > dividend){
            if(skip_flag == 0) skips++;
            quotient<<=1;
        } else {
            skip_flag = 1;
            quotient = (quotient<<1 | 1);
            dividend -= divisor;
        }
    }

    ft_r.characteristic = FLIB_EXTRACT_PART(quotient, FLIB_CARRY_MASK, FLIB_CARRY_SHIFT);
    ft_r.mantissa = FLIB_MASKOFF(quotient, FLIB_MANTISSA_MASK);
    ft_r.exponent = n_exponent- skips;
    ft_r.sign = (ft_a->sign + ft_b->sign) % 2;
    return flib_pack(&ft_r);
}

float __aeabi_i2f(int i) {
    flib_float_t ft_r;

    ft_r.mantissa = 0;
    ft_r.sign = 0;
    ft_r.exponent = 127;

    if( i == 0 ) {
        ft_r.characteristic = 0;
        ft_r.exponent = 0;
    } else if( i > 0 ) {
        ft_r.characteristic = i;
    } else {
        ft_r.characteristic = -i;
        ft_r.sign = 1;
    }

    return flib_pack(&ft_r);
}

float __aeabi_ui2f(unsigned int i) {
    flib_float_t ft_r;

    ft_r.mantissa = 0;
    ft_r.sign = 0;
    ft_r.exponent = 127;

    if( i == 0 ) {
        ft_r.characteristic = 0;
        ft_r.exponent = 0;
    } else {
        ft_r.characteristic = i;
        ft_r.sign = 0;
    }

    return flib_pack(&ft_r);
}

void ftoa(float f, char* buffer) {

    const char _ascii_rep[] = "0123456789";

    // Space to hold the converted mantissa
    char __mantissa_buffer[65];
    char* __mantissa_buffer_ptr = (char*)(__mantissa_buffer);
    __mantissa_buffer_ptr+=64;

    // Space to hold the converted characteristic
    char __characteristic_buffer[33];
    char* __characteristic_buffer_ptr = (char*)(__characteristic_buffer);
    __characteristic_buffer_ptr+=32;

    // Break appart the floating point number and adjust the radix to 0
    // this will put our whole number in characteristic and our rational number in mantissa
    flib_float_t ft_r;
    flib_explode(f, &ft_r);
    flib_adjust_radix(&ft_r, -(ft_r.exponent-127), false);

    uint64_t value = 0, base;

    // Convert the mantissa to a integer representation.
    // The mantissa is  power of 2 fraction. -1**2, -2**2, -3**2, ect.
    for(base = 5000000000000000000ULL; ft_r.mantissa > 0 && base > 0; ft_r.mantissa<<=1, base>>=1) {
        if((ft_r.mantissa & 0x400000) == 0x400000)
            value += base;
    }

    // Convert the digits to ascii. Start at end of char pointer and work backwords.
    *__mantissa_buffer_ptr = '\0';
    __mantissa_buffer_ptr--;

    // To know how many zeros to pad the output with, we assume 19 digits and decrement for each digit we place.
    // The remain value in leading zeros will be the number of zeros that we need to pad.
    int leading_zeros = 19;

    do {
        *__mantissa_buffer_ptr = _ascii_rep[value % 10];
        __mantissa_buffer_ptr--;
        leading_zeros--;
    } while (value /= 10);

    // padd in the leading zeros on the mantissa
    while(leading_zeros-- > 0) {
        *__mantissa_buffer_ptr = '0';
        __mantissa_buffer_ptr--;
    }

    // Convert the characteristic to ascii. Start at end of char pointer and work backwords.
    *__characteristic_buffer_ptr = '\0';
    __characteristic_buffer_ptr--;
    do {
        *__characteristic_buffer_ptr = _ascii_rep[ft_r.characteristic % 10];
        __characteristic_buffer_ptr--;
    } while (ft_r.characteristic /= 10);

    // Setup output buffer. Add sign if needed.
    if(ft_r.sign == 1) { 
        buffer[0] = '-';
        buffer[1] = '\0';
    } else {
        buffer[0] = '\0';
    }

    // Copy in the characteristic and the mantissa
    strncat(buffer, ++__characteristic_buffer_ptr, 32);
    strcat(buffer, ".");
    strncat(buffer, ++__mantissa_buffer_ptr, 6);
}

float atof(char* buffer) {

    flib_float_t ft_r;

    char __buffer[128];
    char *__characteristic_buffer;
    char *__mantissa_buffer;
    int i, offset;

    ft_r.characteristic = 0;
    ft_r.mantissa = 0;
    ft_r.exponent = 127;
    ft_r.sign = 0;

    strncpy(__buffer, buffer, 128);
    __characteristic_buffer = strtok(__buffer, ".");
    __mantissa_buffer = strtok(NULL, ".");

    for(i = strlen(__characteristic_buffer) - 1, offset = 0; i >= 0; i--, offset++) {
        if(__characteristic_buffer[i] == '-') { ft_r.sign = 1; }
        else if(__characteristic_buffer[i] == '+') {}
        else { ft_r.characteristic += (__characteristic_buffer[i] - 0x30) * __uipow(10, offset); }
    }

    float ret = flib_pack(&ft_r); float part = 0;
    for(i = 0; i < strlen(__mantissa_buffer); i++) {

        // Calcualte the value of this digit : (CHAR-0x30) / (10**(i+1))
        part = __aeabi_fdiv(
                __aeabi_i2f(__mantissa_buffer[i] - '0'),
                __aeabi_i2f(__uipow(10, i+1))
        );

        // Add the value to the running total. Invert the addend if the target is negative
        ret = __aeabi_fadd( ret, (ft_r.sign == 1 ? -part : part));
    }

    return ret;
}


#include <math.h>
int main() {

    char ascii_float[] = "123.0001";
    float f = atof(ascii_float);
    float cf = 123.0001;

    printf("Value was:            %s\n"
           "Value interprated as: %f\n"
           "System Sayse:         %f\n",
            ascii_float, f, cf);
    return 0;

    float theF = 44.0217781067;
    char buffer[255] = {0};
    ftoa(theF, buffer);
    printf("Got:      %s\n", buffer);
    printf("Expected: %f\n", theF);
    return 0;

    //float f1 = 16.0;
    //float f2 = 0.00001; //**/
    //float f1 = 200000.5; //**/
    //float f2 = -0.000025;
    //float f1 = 5.0; //**/
    //float f2 = 1.75; //**/
    //float f1 = 1;
    //float f2 = 0;
    float f1 =  __aeabi_ui2f(10u);
    float f2 = -0.1;

    float fa = __aeabi_fdiv(f1, f2);
    float fsa = f1 / f2;
    printf("Divide: %f / %f\n"
           "Got      %f\n"
           "Expected %f\n", f1, f2, fa, fsa);

    print_float(fa);
    print_float(fsa);

    return 0;
}

