INSTALL_PREFIX=/usr/local
BIN_INSTALL_PREFIX=$(INSTALL_PREFIX)/bin

all: install

clean:
	    git clean -f .
install:
	install -m 0755 download_and_format_database.sh $(BIN_INSTALL_PREFIX)
	install -m 0755 create_reads.py $(BIN_INSTALL_PREFIX)
	install -m 0755 insert_SNVs.py $(BIN_INSTALL_PREFIX)
	install -m 0755 insert_SNVs_select_indels.py $(BIN_INSTALL_PREFIX)
	install -m 0755 insert_SVs_and_indels.py $(BIN_INSTALL_PREFIX)
	install -m 0755 insert_SVs.py $(BIN_INSTALL_PREFIX)
	install -m 0755 simulate_SV_BED.py $(BIN_INSTALL_PREFIX)
	install -m 0755 split_fasta_by_contigs.py $(BIN_INSTALL_PREFIX)
.PHONY: install
