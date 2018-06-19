#ifndef MSG_H
#define MSG_H 1

enum LogLevel {msg_debug, msg_verbose, msg_info, msg_warn, msg_error, msg_fatal, msg_silent};

void msg_printf(const enum LogLevel msg_level, const char *fmt, ...);
void msg_abort(const int errret, const char *fmt, ...);

#endif
