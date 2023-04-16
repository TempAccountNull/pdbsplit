// Copyright 2011-2022, Molecular Matters GmbH <office@molecular-matters.com>
// See LICENSE.txt for licensing details (2-clause BSD License: https://opensource.org/licenses/BSD-2-Clause)

#include "PDB_PCH.h"
#include "PDB_RawFile.h"
#include "PDB_Types.h"
#include "PDB_Util.h"
#include "PDB_DirectMSFStream.h"
#include "Foundation/PDB_PointerUtil.h"
#include "Foundation/PDB_Memory.h"
#include "Foundation/PDB_Assert.h"


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
PDB::RawFile::RawFile(RawFile&& other) PDB_NO_EXCEPT
	: m_data(PDB_MOVE(other.m_data))
	, m_blockSize(PDB_MOVE(other.m_blockSize))
	, m_blockCount(PDB_MOVE(other.m_blockCount))
	, m_directoryStream(PDB_MOVE(other.m_directoryStream))
	, m_streamCount(PDB_MOVE(other.m_streamCount))
	, m_streamSizes(PDB_MOVE(other.m_streamSizes))
	, m_streamBlockIndices(PDB_MOVE(other.m_streamBlockIndices))
	, m_streamBlocks(PDB_MOVE(other.m_streamBlocks))
{
	other.m_data = nullptr;
	other.m_blockSize = 0u;
	other.m_blockCount = 0u;
	other.m_streamCount = 0u;
	other.m_streamSizes = nullptr;
	other.m_streamBlockIndices = nullptr;
	other.m_streamBlocks = nullptr;
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
PDB::RawFile& PDB::RawFile::operator=(RawFile&& other) PDB_NO_EXCEPT
{
	if (this != &other)
	{
		PDB_DELETE_ARRAY(m_streamBlocks);

		m_data = PDB_MOVE(other.m_data);
		m_blockSize = PDB_MOVE(other.m_blockSize);
		m_blockCount = PDB_MOVE(other.m_blockCount);
		m_directoryStream = PDB_MOVE(other.m_directoryStream);
		m_streamCount = PDB_MOVE(other.m_streamCount);
		m_streamSizes = PDB_MOVE(other.m_streamSizes);
		m_streamBlockIndices = PDB_MOVE(other.m_streamBlockIndices);
		m_streamBlocks = PDB_MOVE(other.m_streamBlocks);

		other.m_data = nullptr;
		other.m_streamCount = 0u;
		other.m_streamSizes = nullptr;
		other.m_streamBlockIndices = nullptr;
		other.m_streamBlocks = nullptr;
	}

	return *this;
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
PDB::RawFile::RawFile(const void* data) PDB_NO_EXCEPT
	: m_data(data)
	, m_blockSize(0u)
	, m_blockCount(0u)
	, m_directoryStream()
	, m_streamCount(0u)
	, m_streamSizes(nullptr)
	, m_streamBlockIndices(nullptr)
	, m_streamBlocks(nullptr)
{
	bool oldFormat = (std::memcmp(m_data, SuperBlock::MAGIC, sizeof(SuperBlock::MAGIC)) != 0);

	if (oldFormat)
	{
		const SuperBlockV60* superBlock = Pointer::Offset<const SuperBlockV60*>(data, 0u);
		m_blockSize = superBlock->blockSize;
		m_blockCount = superBlock->blockCount;

		const uint32_t directoryBlockCount = PDB::ConvertSizeToBlockCount(superBlock->directorySize, m_blockSize);

		// convert the 16 bit indices
		uint32_t* directoryIndices = PDB_NEW_ARRAY(uint32_t, directoryBlockCount);
		for (uint32_t i = 0u; i < directoryBlockCount; ++i)
		{
			directoryIndices[i] = superBlock->directoryIndices[i];
		}

		m_directoryStream = CoalescedMSFStream(data, m_blockSize, directoryIndices, superBlock->directorySize);
		PDB_DELETE_ARRAY(directoryIndices);
	}
	else
	{
		const SuperBlock* superBlock = Pointer::Offset<const SuperBlock*>(data, 0u);
		m_blockSize = superBlock->blockSize;
		m_blockCount = superBlock->blockCount;

		// the SuperBlock stores an array of indices of blocks that make up the indices of directory blocks, which need to be stitched together to form the directory.
		// the blocks holding the indices of directory blocks are not necessarily contiguous, so they need to be coalesced first.
		const uint32_t directoryBlockCount = PDB::ConvertSizeToBlockCount(superBlock->directorySize, m_blockSize);

		// the directory is made up of directoryBlockCount blocks, so we need that many indices to be read from the blocks that make up the indices
		CoalescedMSFStream directoryIndicesStream(data, m_blockSize, superBlock->directoryBlockIndices, directoryBlockCount * sizeof(uint32_t));

		// these are the indices of blocks making up the directory stream, now guaranteed to be contiguous
		const uint32_t* directoryIndices = directoryIndicesStream.GetDataAtOffset<uint32_t>(0u);

		m_directoryStream = CoalescedMSFStream(data, m_blockSize, directoryIndices, superBlock->directorySize);
	}

	// https://llvm.org/docs/PDB/MsfFile.html#the-stream-directory
	// parse the directory from its contiguous version. the directory matches the following struct:
	//	struct StreamDirectory
	//	{
	//		uint32_t streamCount;
	//		uint32_t streamSizes[streamCount];
	//		uint32_t streamBlocks[streamCount][];
	//	};
	m_streamCount = *m_directoryStream.GetDataAtOffset<uint32_t>(0u);

	m_streamSizes = PDB_NEW_ARRAY(uint32_t, m_streamCount);
	m_streamBlockIndices = PDB_NEW_ARRAY(uint32_t, m_blockCount);

	if (!oldFormat)
	{
		// fast path, just copy the data
		memcpy(m_streamSizes, m_directoryStream.GetDataAtOffset<uint32_t>(sizeof(uint32_t)), sizeof(uint32_t) * m_streamCount);
		memcpy(m_streamBlockIndices, m_directoryStream.GetDataAtOffset<uint32_t>(sizeof(uint32_t) + sizeof(uint32_t) * m_streamCount), sizeof(uint32_t) * m_blockCount);
	}
	else
	{
		// upgrade old format
		struct FileStreamSize
		{
			uint32_t streamSize;
			uint32_t reserved;
		};

		const FileStreamSize* streamSizes = m_directoryStream.GetDataAtOffset<FileStreamSize>(sizeof(uint32_t));
		for (uint32_t i = 0u; i < m_streamCount; ++i)
		{
			m_streamSizes[i] = streamSizes[i].streamSize;
		}

		const uint16_t* streamBlockIndices = m_directoryStream.GetDataAtOffset<uint16_t>(sizeof(uint32_t) + sizeof(FileStreamSize) * m_streamCount);
		for (uint32_t i = 0u; i < m_blockCount; ++i)
		{
			m_streamBlockIndices[i] = streamBlockIndices[i];
		}
	}

	// prepare indices for directly accessing individual streams
	m_streamBlocks = PDB_NEW_ARRAY(const uint32_t*, m_streamCount);

	const uint32_t* indicesForCurrentBlock = m_streamBlockIndices;
	for (uint32_t i = 0u; i < m_streamCount; ++i)
	{
		const uint32_t sizeInBytes = m_streamSizes[i];
		const uint32_t blockCount = ConvertSizeToBlockCount(sizeInBytes, m_blockSize);
		m_streamBlocks[i] = indicesForCurrentBlock;
		
		indicesForCurrentBlock += blockCount;
	}
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
PDB::RawFile::~RawFile(void) PDB_NO_EXCEPT
{
	PDB_DELETE_ARRAY(m_streamBlockIndices);
	PDB_DELETE_ARRAY(m_streamSizes);
	PDB_DELETE_ARRAY(m_streamBlocks);
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
template <typename T>
PDB_NO_DISCARD T PDB::RawFile::CreateMSFStream(uint32_t streamIndex) const PDB_NO_EXCEPT
{
	return T(m_data, m_blockSize, m_streamBlocks[streamIndex], m_streamSizes[streamIndex]);
}


// ------------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------------
template <typename T>
PDB_NO_DISCARD T PDB::RawFile::CreateMSFStream(uint32_t streamIndex, uint32_t streamSize) const PDB_NO_EXCEPT
{
	PDB_ASSERT(streamSize <= m_streamSizes[streamIndex], "Invalid stream size.");

	return T(m_data, m_blockSize, m_streamBlocks[streamIndex], streamSize);
}


// explicit template instantiation
template PDB::CoalescedMSFStream PDB::RawFile::CreateMSFStream<PDB::CoalescedMSFStream>(uint32_t streamIndex) const PDB_NO_EXCEPT;
template PDB::DirectMSFStream PDB::RawFile::CreateMSFStream<PDB::DirectMSFStream>(uint32_t streamIndex) const PDB_NO_EXCEPT;

template PDB::CoalescedMSFStream PDB::RawFile::CreateMSFStream<PDB::CoalescedMSFStream>(uint32_t streamIndex, uint32_t streamSize) const PDB_NO_EXCEPT;
template PDB::DirectMSFStream PDB::RawFile::CreateMSFStream<PDB::DirectMSFStream>(uint32_t streamIndex, uint32_t streamSize) const PDB_NO_EXCEPT;
