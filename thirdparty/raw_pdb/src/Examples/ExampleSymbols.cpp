// Copyright 2011-2022, Molecular Matters GmbH <office@molecular-matters.com>
// See LICENSE.txt for licensing details (2-clause BSD License: https://opensource.org/licenses/BSD-2-Clause)

#include "Examples_PCH.h"
#include "ExampleTimedScope.h"
#include "PDB_RawFile.h"
#include "PDB_DBIStream.h"


namespace
{
	// we don't have to store std::string in the symbols, since all the data is memory-mapped anyway.
	// we do it in this example to ensure that we don't "cheat" when reading the PDB file. memory-mapped data will only
	// be faulted into the process once it's touched, so actually copying the string data makes us touch the needed data,
	// giving us a real performance measurement.
	struct Symbol
	{
		std::string name;
		uint32_t rva;
	};
}


void ExampleSymbols(const PDB::RawFile& rawPdbFile, const PDB::DBIStream& dbiStream);
void ExampleSymbols(const PDB::RawFile& rawPdbFile, const PDB::DBIStream& dbiStream)
{
	TimedScope total("\nRunning example \"Symbols\"");

	// in order to keep the example easy to understand, we load the PDB data serially.
	// note that this can be improved a lot by reading streams concurrently.

	// prepare the image section stream first. it is needed for converting section + offset into an RVA
	TimedScope sectionScope("Reading image section stream");
	const PDB::ImageSectionStream imageSectionStream = dbiStream.CreateImageSectionStream(rawPdbFile);
	sectionScope.Done();


	// prepare the module info stream for matching contributions against files
	TimedScope moduleScope("Reading module info stream");
	const PDB::ModuleInfoStream moduleInfoStream = dbiStream.CreateModuleInfoStream(rawPdbFile);
	moduleScope.Done();


	// prepare symbol record stream needed by both public and global streams
	TimedScope symbolStreamScope("Reading symbol record stream");
	const PDB::CoalescedMSFStream symbolRecordStream = dbiStream.CreateSymbolRecordStream(rawPdbFile);
	symbolStreamScope.Done();

	std::vector<Symbol> symbols;

	// read public symbols
	TimedScope publicScope("Reading public symbol stream");
	const PDB::PublicSymbolStream publicSymbolStream = dbiStream.CreatePublicSymbolStream(rawPdbFile);
	publicScope.Done();
	{
		TimedScope scope("Storing public symbols");

		// incorrect count grabbed here?
		const PDB::ArrayView<PDB::HashRecord> hashRecords = publicSymbolStream.GetRecords();
		//const size_t count = hashRecords.GetLength();

		//symbols.reserve(count);

		for (const PDB::HashRecord& hashRecord : hashRecords)
		{
			const PDB::CodeView::DBI::Record* record = publicSymbolStream.GetRecord(symbolRecordStream, hashRecord);
			if (!record)
			{
				// end of data?
				break;
			}

			const uint32_t rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_PUB32.section, record->data.S_PUB32.offset);
			if (rva == 0u)
			{
				// certain symbols (e.g. control-flow guard symbols) don't have a valid RVA, ignore those
				continue;
			}

			// it's so lame that alloca isn't scope-wide
			char* stringCopy = new char[static_cast<size_t>(record->data.S_PUB32.name.vc60.length) + 1];
			stringCopy[record->data.S_PUB32.name.vc60.length] = '\0';
			memcpy(stringCopy, record->data.S_PUB32.name.vc60.string, record->data.S_PUB32.name.vc60.length);
			symbols.push_back(Symbol { stringCopy, rva });
			delete[] stringCopy;
		}

		//scope.Done(count);
		scope.Done(symbols.size());
	}


	// read global symbols
	TimedScope globalScope("Reading global symbol stream");
	const PDB::GlobalSymbolStream globalSymbolStream = dbiStream.CreateGlobalSymbolStream(rawPdbFile);
	globalScope.Done();
	{
		TimedScope scope("Storing global symbols");

		const PDB::ArrayView<PDB::HashRecord> hashRecords = globalSymbolStream.GetRecords();
		//const size_t count = hashRecords.GetLength();

		//symbols.reserve(symbols.size() + count);
		const size_t old_size = symbols.size();

		for (const PDB::HashRecord& hashRecord : hashRecords)
		{
			const PDB::CodeView::DBI::Record* record = globalSymbolStream.GetRecord(symbolRecordStream, hashRecord);
			if (!record)
			{
				// end of data?
				break;
			}

			const char* name = nullptr;
			uint32_t rva = 0u;
			if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_GDATA32)
			{
				name = record->data.S_GDATA32.name.vc70.string;
				rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_GDATA32.section, record->data.S_GDATA32.offset);
			}
			else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_GTHREAD32)
			{
				name = record->data.S_GTHREAD32.name.vc70.string;
				rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_GTHREAD32.section, record->data.S_GTHREAD32.offset);
			}
			else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LDATA32)
			{
				name = record->data.S_LDATA32.name.vc70.string;
				rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LDATA32.section, record->data.S_LDATA32.offset);
			}
			else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LTHREAD32)
			{
				name = record->data.S_LTHREAD32.name.vc70.string;
				rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LTHREAD32.section, record->data.S_LTHREAD32.offset);
			}
			else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_UDT)
			{
				name = record->data.S_UDT.name.vc70.string;
			}
			else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_UDT_ST)
			{
				name = record->data.S_UDT_ST.name.vc70.string;
			}

			if (rva == 0u)
			{
				// certain symbols (e.g. control-flow guard symbols) don't have a valid RVA, ignore those
				continue;
			}

			symbols.push_back(Symbol { name, rva });
		}

		scope.Done(symbols.size() - old_size);
	}


	// read module symbols
	{
		TimedScope scope("Storing symbols from modules");

		const PDB::ArrayView<PDB::ModuleInfoStream::Module> modules = moduleInfoStream.GetModules();

		for (const PDB::ModuleInfoStream::Module& module : modules)
		{
			if (!module.HasSymbolStream())
			{
				continue;
			}

			const PDB::ModuleSymbolStream moduleSymbolStream = module.CreateSymbolStream(rawPdbFile);
			moduleSymbolStream.ForEachSymbol([&symbols, &imageSectionStream](const PDB::CodeView::DBI::Record* record)
			{
				const char* name = nullptr;
				uint32_t rva = 0u;
				if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_THUNK32)
				{
					if (record->data.S_THUNK32.thunk == PDB::CodeView::DBI::ThunkOrdinal::TrampolineIncremental)
					{
						// we have never seen incremental linking thunks stored inside a S_THUNK32 symbol, but better be safe than sorry
						name = "ILT";
						rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_THUNK32.section, record->data.S_THUNK32.offset);
					}
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_TRAMPOLINE)
				{
					// incremental linking thunks are stored in the linker module
					name = "ILT";
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_TRAMPOLINE.thunkSection, record->data.S_TRAMPOLINE.thunkOffset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_BLOCK32)
				{
					// blocks never store a name and are only stored for indicating whether other symbols are children of this block
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LABEL32)
				{
					// labels don't have a name
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LPROC32)
				{
					name = record->data.S_LPROC32.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LPROC32.section, record->data.S_LPROC32.offset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_GPROC32)
				{
					name = record->data.S_GPROC32.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_GPROC32.section, record->data.S_GPROC32.offset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LPROC32_ID)
				{
					name = record->data.S_LPROC32_ID.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LPROC32_ID.section, record->data.S_LPROC32_ID.offset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_GPROC32_ID)
				{
					name = record->data.S_GPROC32_ID.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_GPROC32_ID.section, record->data.S_GPROC32_ID.offset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LDATA32)
				{
					name = record->data.S_LDATA32.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LDATA32.section, record->data.S_LDATA32.offset);
				}
				else if (record->header.kind == PDB::CodeView::DBI::SymbolRecordKind::S_LTHREAD32)
				{
					name = record->data.S_LTHREAD32.name.vc70.string;
					rva = imageSectionStream.ConvertSectionOffsetToRVA(record->data.S_LTHREAD32.section, record->data.S_LTHREAD32.offset);
				}

				if (rva == 0u)
				{
					// certain symbols (e.g. control-flow guard symbols) don't have a valid RVA, ignore those
					return;
				}

				symbols.push_back(Symbol { name, rva });
			});
		}

		scope.Done(modules.GetLength());
	}

	TimedScope sortScope("std::sort symbols");
	std::sort(symbols.begin(), symbols.end(), [](const Symbol& lhs, const Symbol& rhs)
		{
			return lhs.rva < rhs.rva;
		});
	sortScope.Done();

	total.Done(symbols.size());

	printf("Symbol list:\n");
	for (const Symbol& symbol : symbols)
	{
		printf("%08x : %s\n", symbol.rva, symbol.name.c_str());
	}
}
