
class O2scl < Formula
  homepage "http://o2scl.sourceforge.net"
  url "http://23.115.64.197/o2scl-0.918.tar.gz"
  sha256 "2b5ee94f937472a1a99103f6319f3bebb9f9e0f9118d2870e5f8177ac999586a"

  depends_on "gsl"
  depends_on "hdf5"
  depends_on "boost"
  depends_on "readline"

  def install
    system "./configure", "--disable-debug",
                          "--disable-dependency-tracking",
                          "--disable-silent-rules",
                          "--prefix=#{prefix}"
    system "make"
    system "make", "install"
  end

  test do
    system "acol --version"
  end
end
