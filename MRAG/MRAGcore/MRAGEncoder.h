#ifndef log2
#define log2(x) (log((double)(x)) / (double)log(2.))
#endif
namespace MRAG {

template <typename DataType> class Encoder {
protected:
  static const bool bVerbose = false;
  bool m_bEncoded;
  unsigned int m_nQuantization;
  unsigned int m_nItems;
  BitStream m_encodedData;

public:
  Encoder()
      : m_nQuantization(0), m_encodedData(), m_bEncoded(false), m_nItems(0) {}

  ~Encoder() { m_encodedData.dispose(); }

  virtual void encode(const std::vector<DataType> &stream, const int nSymbols) {
    assert(m_bEncoded == false);
    assert(nSymbols >= 2);

    m_nQuantization =
        std::max((unsigned int)0, (unsigned int)ceil(log2(nSymbols)));
    m_nItems = stream.size();

    const double uncompressedMB = m_nItems * sizeof(DataType) / 1024. / 1024.;
    const double compressionFactor = uncompressedMB / getMemorySize();

    if (bVerbose)
      printf("Encoder Compression Rate: %f (%.2f KB instead of %.2f KB)\n",
             compressionFactor, getMemorySize() * 1024, uncompressedMB * 1024.);

    m_encodedData.setup(m_nItems * m_nQuantization);

    const typename std::vector<DataType>::const_iterator itE = stream.end();
    for (typename std::vector<DataType>::const_iterator itSource =
             stream.begin();
         itSource != itE; itSource++)
      m_encodedData.append_bits((unsigned int)*itSource, m_nQuantization);

    m_bEncoded = true;
  }

  virtual void decode(std::vector<DataType> &stream) const {
    assert(m_bEncoded == true);
    assert(stream.size() == 0);

    stream.resize(m_nItems);
    const int n = m_nItems;
    typename std::vector<DataType>::iterator itDest = stream.begin();
    for (int i = 0; i < n; i++, itDest++)
      *itDest = (DataType)m_encodedData.get_bits(m_nQuantization * i,
                                                 m_nQuantization);
  }

  virtual float getMemorySize() const {
    return (sizeof(Encoder) + m_encodedData.getLength() / 8.) / (1024. * 1024.);
  }

private:
  // forbidden
  Encoder(const Encoder &h)
      : m_nQuantization(0), m_encodedData(), m_bEncoded(false), m_nItems(0) {
    abort();
  }
  Encoder &operator=(const Encoder &h) {
    abort();
    return *this;
  }
};

} // namespace MRAG
