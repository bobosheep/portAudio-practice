main:	main.c
	gcc main.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o main
pa_devs: 	pa_devs.c
	gcc pa_devs.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o pa_devs
paex_read_write_wire: 	paex_read_write_wire.c
	gcc paex_read_write_wire.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o paex_read_write_wire
paex_record: 	paex_record.c
	gcc paex_record.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o paex_record
myRecord: myRecord.c
	gcc myRecord.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o myRecord
myRecord2: myRecord2.c
	gcc myRecord2.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o myRecord2
myStreamtest: myStreamtest.c
	gcc myStreamtest.c libportaudio.a -lrt -lm -lasound -ljack -pthread -o myStreamtest