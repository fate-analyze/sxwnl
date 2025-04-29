import math
from algo_lib.DateTime import DateTime
from astropy.coordinates import get_sun, get_body, EarthLocation, CIRS, GeocentricTrueEcliptic
import astropy.units as units
from astropy.coordinates import solar_system_ephemeris
from astropy.time import Time
from datetime import timezone


class TimeUtil:
    @staticmethod
    def datetime_2_julian(dt: DateTime) -> float:
        """
        日期转为儒略日
        :param dt: 时间
        """
        datetime_str = dt.to_str(True)
        t = Time(datetime_str, format='fits', scale='utc')
        return t.jd1 + t.jd2

    @staticmethod
    def julian_2_datetime(jd: float) -> DateTime:
        """
        儒略日转日期
        :param jd: 儒略日
        :return: 日期
        """
        t = Time(jd, format='jd', scale='utc')
        frac, integer = math.modf(t.ymdhms.second)
        dt = DateTime(t.ymdhms.year, t.ymdhms.month, t.ymdhms.day, t.ymdhms.hour, t.ymdhms.minute, int(integer),
                      int(round(frac * 1000)), timezone.utc)

        return dt

    @staticmethod
    def datetime_2_sun(dt : DateTime, lng : float) -> float:
        """
        计算给定时间和经度的太阳视黄经和本地真太阳时
        :param dt: 时间
        :param lng: 经度(单位：度）
        """
        with solar_system_ephemeris.set('jpl'):
            t = Time(dt.to_str(True), format='fits', scale='utc')

            # 计算太阳视黄经
            sun = get_sun(t)
            ecliptic = GeocentricTrueEcliptic(equinox=t)
            sun_ecl = sun.transform_to(ecliptic)
            solar_lon = sun_ecl.lon.degrees  # 太阳视黄经（单位：度）

            # 计算真太阳时
            # 转换到CIRS坐标系获取视赤经
            sun_cirs = sun.transform_to(CIRS(obstime=t))
            ra = sun_cirs.ra

            # 计算本地视恒星时
            lst = t.sidereal_time('apparent', longitude=lng * units.si.deg)

            # 计算时角并转换为小时
            ha = (lst - ra).to(units.si.hourangle)  # 时角（小时）
            true_solar = (ha.value + 12) % 24

            # 四舍五入到最近小时并确保在0-23范围内
            true_solar_hour = int(round(true_solar)) % 24

            return round(solar_lon, 2), true_solar_hour