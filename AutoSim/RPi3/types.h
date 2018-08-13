#pragma once

#if defined(_MSC_VER)
#define _GLIBCXX_ABI_TAG_CXX11
#define __size_t__
#if defined(__cplusplus)
#undef __cplusplus
#define __cplusplus 201402L
#endif

#ifndef __builtin_expect
#define __builtin_expect(x,y) (x)
#endif
#ifndef __ASSERT_FUNCTION
#define __ASSERT_FUNCTION
#endif
#ifndef __builtin_memset
#define __builtin_memset(x,y)
#endif
#ifndef uint8_t
typedef unsigned __int8 uint8_t;
#endif //uint8_t
#ifndef uint16_t
typedef unsigned __int16 uint16_t;
#endif //uint16_t
#ifndef uint32_t
typedef unsigned __int32 uint32_t;
#endif //uint32_t
#ifndef uint64_t
typedef unsigned __int64 uint64_t;
#endif //uint64_t
#endif